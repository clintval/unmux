//! The check of a `--reference` FASTA against each CRAM's `@SQ` `LN` and `M5`.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use anyhow::{anyhow, bail, Context, Result};
use md5::{Digest, Md5};
use noodles::fasta::{self, fai};
use noodles::sam;
use noodles::sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;
use rayon::prelude::*;

/// An MD5 digest.
type Md5Digest = [u8; 16];

/// Small reads contend in the kernel when many contigs are hashed at once.
const READ_BUFFER: usize = 1 << 20;

/// Enough to read the bytes around one end of a contig.
const BOUNDS_BUFFER: usize = 1 << 12;

/// A reader of the `--reference` FASTA, uncompressed or bgzipped.
type Reader = fasta::io::BufReader<File>;

/// Checks CRAM headers against an indexed FASTA, hashing each contig at most once on `threads` threads.
pub(crate) struct HeaderCheck {
    path: PathBuf,
    index: fai::Index,
    digests: Vec<OnceLock<ContigDigest>>,
    threads: usize,
}

impl HeaderCheck {
    /// A check against the FASTA at `path`, whose `.fai` is `index`.
    pub(crate) fn new(path: &Path, index: fai::Index, threads: usize) -> Self {
        Self {
            path: path.to_path_buf(),
            digests: index.as_ref().iter().map(|_| OnceLock::new()).collect(),
            index,
            threads: threads.max(1),
        }
    }

    /// Fail unless every `@SQ` contig is in this FASTA with the same `LN` and any `M5`.
    pub(crate) fn check(&self, header: &sam::Header, input: &Path) -> Result<()> {
        let input = &input_name(input);
        let expected = expected_contigs(header, input)?;
        let positions = self.positions_for(&expected, input)?;
        let (hashed, unhashed): (Vec<_>, Vec<_>) = expected
            .iter()
            .zip(positions)
            .partition(|(contig, _)| contig.md5.is_some());
        if !unhashed.is_empty() {
            log::warn!(
                "{} of {} @SQ contig(s) in {} have no M5, so only their lengths were checked against the --reference FASTA",
                unhashed.len(),
                expected.len(),
                input
            );
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .context("failed to build the reference check thread pool")?;
        let (bounds, digests): (Vec<Result<()>>, Vec<Result<ContigDigest>>) = pool.install(|| {
            let records = self.index.as_ref();
            let bounds = unhashed
                .par_iter()
                .map_init(
                    || None,
                    |reader, &(_, i)| self.check_bounds(reader, &records[i]),
                )
                .collect();
            let digests = hashed
                .par_iter()
                .with_max_len(1)
                .map(|&(_, i)| self.cached_digest(i))
                .collect();
            (bounds, digests)
        });
        let (mut bounds, mut digests) = (bounds.into_iter(), digests.into_iter());
        for contig in &expected {
            let Some(md5) = contig.md5 else {
                bounds
                    .next()
                    .expect("one bounds check per contig without M5")?;
                continue;
            };
            let actual = digests.next().expect("one digest per contig with M5")?;
            if !actual.matches(&md5) {
                bail!(
                    "contig `{}` has MD5 {} in the --reference FASTA {} but @SQ M5:{} in {}; the FASTA is not the one the CRAM was written against",
                    contig.name,
                    to_hex(&actual.spec),
                    self.path.display(),
                    to_hex(&md5),
                    input
                );
            }
        }
        Ok(())
    }

    /// The digest of the `i`th `.fai` record, hashing it only on first use.
    fn cached_digest(&self, i: usize) -> Result<ContigDigest> {
        if let Some(digest) = self.digests[i].get() {
            return Ok(*digest);
        }
        let digest = self.digest(&self.index.as_ref()[i])?;
        Ok(*self.digests[i].get_or_init(|| digest))
    }

    /// The `.fai` position of each contig, failing on one missing or of another length.
    fn positions_for(&self, expected: &[ExpectedContig], input: &str) -> Result<Vec<usize>> {
        let by_name: HashMap<&[u8], (usize, &fai::Record)> = self
            .index
            .as_ref()
            .iter()
            .enumerate()
            .map(|(i, record)| (record.name().as_ref(), (i, record)))
            .collect();
        expected
            .iter()
            .map(|contig| {
                let &(i, record) = by_name.get(contig.name.as_bytes()).ok_or_else(|| {
                    anyhow!(
                        "contig `{}` (@SQ in {}) is missing from the --reference FASTA {}; pass the FASTA the CRAM was written against, or --skip-reference-check if no record is on that contig",
                        contig.name,
                        input,
                        self.path.display()
                    )
                })?;
                if record.length() != contig.length {
                    bail!(
                        "contig `{}` is {} bp in the --reference FASTA {} but @SQ LN:{} in {}; the FASTA is not the one the CRAM was written against",
                        contig.name,
                        record.length(),
                        self.path.display(),
                        contig.length,
                        input
                    );
                }
                Ok(i)
            })
            .collect()
    }

    /// The MD5s of one sequence with bytes outside `!`..=`~` dropped, as [`ContigDigest`] normalizes it.
    fn digest(&self, record: &fai::Record) -> Result<ContigDigest> {
        let mut reader = self.open(READ_BUFFER)?;
        self.seek_to_start(&mut reader, record)?;
        let length = record.length();
        let mut remaining = length;
        let mut hasher = ContigHasher::default();
        let mut normalized = Vec::new();
        while remaining > 0 {
            let buf = reader.fill_buf().with_context(|| self.read_context())?;
            if buf.is_empty() {
                return Err(self.stale(
                    record,
                    format!("has {} bp, not its indexed {length}", length - remaining),
                ));
            }
            normalized.clear();
            let mut consumed = 0;
            while consumed < buf.len() && remaining > 0 {
                let rest = &buf[consumed..];
                let run = rest
                    .iter()
                    .position(|&byte| !byte.is_ascii_graphic() || byte == b'>')
                    .unwrap_or(rest.len());
                let take = run.min(usize::try_from(remaining).unwrap_or(usize::MAX));
                normalized.extend_from_slice(&rest[..take]);
                remaining -= take as u64;
                consumed += take;
                if take == run && run < rest.len() {
                    if rest[run] != b'>' {
                        consumed += 1;
                    } else if remaining > 0 {
                        return Err(self.stale(
                            record,
                            format!("has {} bp, not its indexed {length}", length - remaining),
                        ));
                    }
                }
            }
            hasher.update(&normalized);
            reader.consume(consumed);
        }
        self.expect_end(&mut reader, record)?;
        Ok(hasher.finalize())
    }

    /// Fail unless `record` starts and ends where the `.fai` says, reading only the bytes around its ends.
    fn check_bounds(&self, reader: &mut Option<Reader>, record: &fai::Record) -> Result<()> {
        let reader = match reader {
            Some(reader) => reader,
            None => reader.insert(self.open(BOUNDS_BUFFER)?),
        };
        self.seek_to_start(reader, record)?;
        if let Some(last) = record.length().checked_sub(1) {
            let line_bases = record.line_bases().max(1);
            let position =
                record.offset() + last / line_bases * record.line_width() + last % line_bases;
            reader
                .seek(SeekFrom::Start(position))
                .with_context(|| self.read_context())?;
            let byte = reader
                .fill_buf()
                .with_context(|| self.read_context())?
                .first();
            if !byte.is_some_and(|&byte| byte.is_ascii_graphic() && byte != b'>') {
                return Err(self.stale(
                    record,
                    format!("is shorter than its indexed {} bp", record.length()),
                ));
            }
            reader.consume(1);
        }
        self.expect_end(reader, record)
    }

    /// Open the FASTA for reading, with a read buffer of `capacity` bytes when uncompressed.
    fn open(&self, capacity: usize) -> Result<Reader> {
        let context = || format!("failed to open {}", self.path.display());
        match self.path.extension().and_then(|ext| ext.to_str()) {
            Some("gz" | "bgz") => Ok(fasta::io::indexed_reader::Builder::default()
                .set_index(fai::Index::default())
                .build_from_path(&self.path)
                .with_context(context)?
                .into_inner()),
            _ => Ok(fasta::io::BufReader::Uncompressed(
                BufReader::with_capacity(capacity, File::open(&self.path).with_context(context)?),
            )),
        }
    }

    /// Seek to the first base of `record`, failing unless a line starts there.
    fn seek_to_start(&self, reader: &mut Reader, record: &fai::Record) -> Result<()> {
        let starts_a_line = match record.offset().checked_sub(1) {
            Some(before) => {
                reader
                    .seek(SeekFrom::Start(before))
                    .with_context(|| self.read_context())?;
                reader
                    .fill_buf()
                    .with_context(|| self.read_context())?
                    .first()
                    == Some(&b'\n')
            }
            None => false,
        };
        if !starts_a_line {
            return Err(self.stale(
                record,
                format!("does not start at its indexed offset {}", record.offset()),
            ));
        }
        reader.consume(1);
        Ok(())
    }

    /// Fail unless only line breaks lie between the reader and the next `>` or the end of the file.
    fn expect_end(&self, reader: &mut Reader, record: &fai::Record) -> Result<()> {
        loop {
            let buf = reader.fill_buf().with_context(|| self.read_context())?;
            match buf.iter().position(u8::is_ascii_graphic) {
                Some(i) if buf[i] == b'>' => return Ok(()),
                Some(_) => {
                    return Err(self.stale(
                        record,
                        format!("is longer than its indexed {} bp", record.length()),
                    ))
                }
                None if buf.is_empty() => return Ok(()),
                None => {
                    let skipped = buf.len();
                    reader.consume(skipped);
                }
            }
        }
    }

    /// The error for a `.fai` record that no longer matches the FASTA.
    fn stale(&self, record: &fai::Record, problem: String) -> anyhow::Error {
        let mut index_path = self.path.clone().into_os_string();
        index_path.push(".fai");
        anyhow!(
            "the .fai index {} is stale: contig `{}` {problem} in {}; rebuild it with samtools faidx",
            Path::new(&index_path).display(),
            record.name(),
            self.path.display()
        )
    }

    /// The context of a failed read of the FASTA.
    fn read_context(&self) -> String {
        format!("failed to read {}", self.path.display())
    }
}

/// The MD5s an `@SQ` `M5` may hold: letters uppercased as in the SAM spec, or bit 0x20 cleared as htslib writes CRAM.
#[derive(Clone, Copy, Debug, PartialEq)]
struct ContigDigest {
    spec: Md5Digest,
    htslib: Option<Md5Digest>,
}

impl ContigDigest {
    /// Whether `md5` is either digest.
    fn matches(&self, md5: &Md5Digest) -> bool {
        &self.spec == md5 || self.htslib.as_ref() == Some(md5)
    }
}

/// Hashes a sequence both ways, starting the htslib digest only at its first non-letter byte.
#[derive(Default)]
struct ContigHasher {
    spec: Md5,
    htslib: Option<Md5>,
    buf: Vec<u8>,
}

impl ContigHasher {
    /// Add the graphic bytes `bases` to both digests.
    fn update(&mut self, mut bases: &[u8]) {
        if self.htslib.is_none() {
            let letters = bases
                .iter()
                .position(|byte| !byte.is_ascii_alphabetic())
                .unwrap_or(bases.len());
            self.update_spec(&bases[..letters]);
            if letters == bases.len() {
                return;
            }
            self.htslib = Some(self.spec.clone());
            bases = &bases[letters..];
        }
        self.update_spec(bases);
        if let Some(htslib) = &mut self.htslib {
            self.buf.clear();
            self.buf.extend(bases.iter().map(|byte| byte & !0x20));
            htslib.update(&self.buf);
        }
    }

    /// Add `bases` to the SAM spec digest.
    fn update_spec(&mut self, bases: &[u8]) {
        self.buf.clear();
        self.buf.extend_from_slice(bases);
        self.buf.make_ascii_uppercase();
        self.spec.update(&self.buf);
    }

    /// Both digests, the htslib one only when it can differ.
    fn finalize(self) -> ContigDigest {
        ContigDigest {
            spec: self.spec.finalize().into(),
            htslib: self.htslib.map(|htslib| htslib.finalize().into()),
        }
    }
}

/// One `@SQ` line's name, length and `M5`.
struct ExpectedContig {
    name: String,
    length: u64,
    md5: Option<Md5Digest>,
}

/// The `@SQ` lines of `header`, in header order.
fn expected_contigs(header: &sam::Header, input: &str) -> Result<Vec<ExpectedContig>> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, map)| {
            let name = name.to_string();
            let md5 = map
                .other_fields()
                .get(&MD5_CHECKSUM)
                .map(|value| {
                    parse_md5(value).ok_or_else(|| {
                        anyhow!(
                            "contig `{name}` has a malformed @SQ M5:{value} in {} (expected 32 hex digits)",
                            input
                        )
                    })
                })
                .transpose()?;
            Ok(ExpectedContig {
                length: map.length().get() as u64,
                md5,
                name,
            })
        })
        .collect()
}

/// An input path for messages, with `-` read as stdin.
fn input_name(input: &Path) -> String {
    if input.as_os_str() == "-" {
        "stdin".to_string()
    } else {
        input.display().to_string()
    }
}

/// Parse a 32-digit hex MD5, in either case.
fn parse_md5(hex: &[u8]) -> Option<Md5Digest> {
    if hex.len() != 32 || !hex.iter().all(u8::is_ascii_hexdigit) {
        return None;
    }
    let hex = std::str::from_utf8(hex).ok()?;
    u128::from_str_radix(hex, 16).ok().map(u128::to_be_bytes)
}

/// Lowercase hex, as `@SQ M5` writes it.
fn to_hex(digest: &Md5Digest) -> String {
    digest.iter().map(|byte| format!("{byte:02x}")).collect()
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use noodles::sam::header::record::value::map::ReferenceSequence;
    use noodles::sam::header::record::value::Map;

    use super::*;

    const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/ref.fa");
    const FIXTURE_MD5: &str = "e1a82f39b96c07c57c14cfb6ec9e0caa";

    fn header(contigs: &[(&str, usize, Option<&str>)]) -> sam::Header {
        let mut builder = sam::Header::builder();
        for &(name, length, md5) in contigs {
            let mut map = Map::<ReferenceSequence>::new(NonZero::new(length).unwrap());
            if let Some(md5) = md5 {
                map.other_fields_mut().insert(MD5_CHECKSUM, md5.into());
            }
            builder = builder.add_reference_sequence(name, map);
        }
        builder.build()
    }

    fn open(path: &Path, threads: usize) -> HeaderCheck {
        let index = fai::fs::read(format!("{}.fai", path.display())).unwrap();
        HeaderCheck::new(path, index, threads)
    }

    fn check(contigs: &[(&str, usize, Option<&str>)]) -> Result<()> {
        open(Path::new(FIXTURE), 2).check(&header(contigs), Path::new("in.cram"))
    }

    fn write_fasta(dir: &Path, body: &str) -> PathBuf {
        let path = dir.join("ref.fa");
        std::fs::write(&path, body).unwrap();
        fai::fs::write(dir.join("ref.fa.fai"), &fasta::fs::index(&path).unwrap()).unwrap();
        path
    }

    #[test]
    fn test_matching_reference_passes() {
        check(&[("chr1", 3000, Some(FIXTURE_MD5))]).unwrap();
        check(&[("chr1", 3000, Some(&FIXTURE_MD5.to_uppercase()))]).unwrap();
    }

    #[test]
    fn test_contig_without_m5_checks_length_only() {
        check(&[("chr1", 3000, None)]).unwrap();
        let err = check(&[("chr1", 2999, None)]).unwrap_err().to_string();
        assert!(err.contains("contig `chr1` is 3000 bp"), "{err}");
        assert!(err.contains("@SQ LN:2999 in in.cram"), "{err}");
    }

    #[test]
    fn test_md5_mismatch_fails() {
        let wrong = "0".repeat(32);
        let err = check(&[("chr1", 3000, Some(&wrong))])
            .unwrap_err()
            .to_string();
        assert!(err.contains("contig `chr1` has MD5 e1a82f39"), "{err}");
        assert!(err.contains(&format!("@SQ M5:{wrong} in in.cram")), "{err}");
    }

    #[test]
    fn test_missing_contig_fails() {
        let err = check(&[("chr1", 3000, None), ("chrM", 16569, None)])
            .unwrap_err()
            .to_string();
        assert!(err.contains("contig `chrM`"), "{err}");
        assert!(err.contains("missing from the --reference FASTA"), "{err}");
        assert!(err.contains("--skip-reference-check"), "{err}");
    }

    #[test]
    fn test_malformed_m5_fails() {
        let err = check(&[("chr1", 3000, Some("e1a8"))])
            .unwrap_err()
            .to_string();
        assert!(err.contains("malformed @SQ M5:e1a8"), "{err}");
    }

    #[test]
    fn test_stdin_input_is_named_stdin() {
        let err = open(Path::new(FIXTURE), 1)
            .check(&header(&[("chrM", 16569, None)]), Path::new("-"))
            .unwrap_err()
            .to_string();
        assert!(err.contains("(@SQ in stdin)"), "{err}");
    }

    #[test]
    fn test_header_without_sq_passes() {
        check(&[]).unwrap();
    }

    #[test]
    fn test_digest_matches_htslib_normalization() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">a\nacgT\r\nNN\r\n>b desc\nAC\n");
        let check = open(&path, 1);
        let records = check.index.as_ref();
        let digest = |i: usize| to_hex(&check.digest(&records[i]).unwrap().spec);
        assert_eq!(digest(0), to_hex(&Md5::digest(b"ACGTNN").into()));
        assert_eq!(digest(1), to_hex(&Md5::digest(b"AC").into()));
    }

    #[test]
    fn test_m5_from_htslib_cram_writer_also_matches() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), &format!(">a\n{}", "ACGT*-.ACG\n".repeat(5)));
        let check = open(&path, 1);
        let samtools_dict = "e6b475a6ea3c95c302e8879b292cb6a7";
        let samtools_view = "c345582aa180e8a916ec7eac268d3e64";
        for md5 in [samtools_dict, samtools_view] {
            check
                .check(&header(&[("a", 50, Some(md5))]), Path::new("in.cram"))
                .unwrap();
        }
        let err = check
            .check(
                &header(&[("a", 50, Some(&"0".repeat(32)))]),
                Path::new("in.cram"),
            )
            .unwrap_err()
            .to_string();
        assert!(err.contains(&format!("has MD5 {samtools_dict}")), "{err}");
    }

    #[test]
    fn test_htslib_digest_is_computed_only_for_non_letter_bytes() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">a\nacgTN\n>b\nAC\nG*\n");
        let check = open(&path, 1);
        let records = check.index.as_ref();
        assert!(check.digest(&records[0]).unwrap().htslib.is_none());
        let digest = check.digest(&records[1]).unwrap();
        assert_eq!(digest.spec, <[u8; 16]>::from(Md5::digest(b"ACG*")));
        assert_eq!(digest.htslib, Some(Md5::digest(b"ACG\n").into()));
    }

    #[test]
    fn test_contig_hasher_is_independent_of_chunking() {
        let bases = b"acgTN*-.09{}acgt";
        let whole = {
            let mut hasher = ContigHasher::default();
            hasher.update(bases);
            hasher.finalize()
        };
        for split in 0..=bases.len() {
            let mut hasher = ContigHasher::default();
            hasher.update(&bases[..split]);
            hasher.update(&bases[split..]);
            assert_eq!(hasher.finalize(), whole, "split at {split}");
        }
        let htslib: Vec<u8> = bases.iter().map(|byte| byte & !0x20).collect();
        assert_eq!(whole.htslib, Some(Md5::digest(&htslib).into()));
    }

    /// Index `before`, then overwrite the FASTA with `after` and digest contig `i` against the old index.
    fn digest_with_stale_index(before: &str, after: &str, i: usize) -> Result<ContigDigest> {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), before);
        std::fs::write(&path, after).unwrap();
        let check = open(&path, 1);
        check.digest(&check.index.as_ref()[i])
    }

    #[test]
    fn test_digest_detects_a_stale_index() {
        let two = ">a\nACGT\n>b\nTT\n";
        let cases = [
            (
                ">a\nACGT\n",
                ">a\nAC",
                0,
                "contig `a` has 2 bp, not its indexed 4",
            ),
            (
                two,
                ">a\nAC\n>b\nTT\n",
                0,
                "contig `a` has 2 bp, not its indexed 4",
            ),
            (
                two,
                ">a\nACGTA\n>b\nTT\n",
                0,
                "contig `a` is longer than its indexed 4 bp",
            ),
            (
                two,
                ">a\nACGT\n>bb\nTT\n",
                1,
                "contig `b` does not start at its indexed offset 11",
            ),
            (
                two,
                ">aa\nACGT\n>b\nTT\n",
                0,
                "contig `a` does not start at its indexed offset 3",
            ),
        ];
        for (before, after, i, expected) in cases {
            let err = digest_with_stale_index(before, after, i)
                .unwrap_err()
                .to_string();
            assert!(err.contains("the .fai index "), "{err}");
            assert!(err.contains("ref.fa.fai is stale"), "{err}");
            assert!(err.contains(expected), "{after:?}: {err}");
        }
    }

    /// Index `before`, then overwrite the FASTA with `after` and check the bounds of contig `i`.
    fn bounds_with_stale_index(before: &str, after: &str, i: usize) -> Result<()> {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), before);
        std::fs::write(&path, after).unwrap();
        let check = open(&path, 1);
        check.check_bounds(&mut None, &check.index.as_ref()[i])
    }

    #[test]
    fn test_bounds_check_detects_a_stale_index() {
        let two = ">a\nACGT\n>b\nTT\n";
        let cases = [
            (
                ">a\nACGT\n",
                ">a\nAC",
                0,
                "contig `a` is shorter than its indexed 4 bp",
            ),
            (
                two,
                ">a\nAC\n>b\nTT\n",
                0,
                "contig `a` is shorter than its indexed 4 bp",
            ),
            (
                two,
                ">a\nACGTA\n>b\nTT\n",
                0,
                "contig `a` is longer than its indexed 4 bp",
            ),
            (
                two,
                ">a\nACGT\n>bb\nTT\n",
                1,
                "contig `b` does not start at its indexed offset 11",
            ),
            (
                two,
                ">aa\nACGT\n>b\nTT\n",
                0,
                "contig `a` does not start at its indexed offset 3",
            ),
        ];
        for (before, after, i, expected) in cases {
            let err = bounds_with_stale_index(before, after, i)
                .unwrap_err()
                .to_string();
            assert!(err.contains("ref.fa.fai is stale"), "{err}");
            assert!(err.contains(expected), "{after:?}: {err}");
        }
    }

    #[test]
    fn test_bounds_check_accepts_wrapped_lines_line_endings_and_eof() {
        for fasta in [
            ">a\nACG\nTAC\nGT\n>b\nTT",
            ">a\r\nACG\r\nTAC\r\nGT\r\n>b\nTT\n\n",
        ] {
            for i in 0..2 {
                bounds_with_stale_index(fasta, fasta, i).unwrap();
            }
        }
    }

    #[test]
    fn test_contig_without_m5_fails_on_a_stale_index() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">a\nACGT\n>b\nTT\n");
        std::fs::write(&path, ">a description\nACGT\n>b\nTT\n").unwrap();
        let md5 = to_hex(&Md5::digest(b"TT").into());
        let err = open(&path, 2)
            .check(
                &header(&[("a", 4, None), ("b", 2, Some(&md5))]),
                Path::new("in.cram"),
            )
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("contig `a` does not start at its indexed offset 3"),
            "{err}"
        );
    }

    #[test]
    fn test_digest_accepts_line_endings_trailing_blank_lines_and_eof() {
        let fasta = ">a\nACGT\r\n>b\nTT";
        let digest = |i| to_hex(&digest_with_stale_index(fasta, fasta, i).unwrap().spec);
        assert_eq!(digest(0), to_hex(&Md5::digest(b"ACGT").into()));
        assert_eq!(digest(1), to_hex(&Md5::digest(b"TT").into()));
        let blank = digest_with_stale_index(">a\nACGT\n", ">a\nACGT\n\n\n", 0).unwrap();
        assert_eq!(to_hex(&blank.spec), to_hex(&Md5::digest(b"ACGT").into()));
    }

    #[test]
    fn test_stale_index_errors_name_the_first_contig_in_header_order() {
        let dir = tempfile::tempdir().unwrap();
        let body: String = (0..8).map(|i| format!(">c{i}\nACGTACGT\n")).collect();
        let path = write_fasta(dir.path(), &body);
        std::fs::write(&path, body.replacen("ACGTACGT", "ACGT", 1)).unwrap();
        let md5 = to_hex(&Md5::digest(b"ACGTACGT").into());
        let names: Vec<String> = (0..8).map(|i| format!("c{i}")).collect();
        let contigs: Vec<(&str, usize, Option<&str>)> = names
            .iter()
            .map(|name| (name.as_str(), 8, Some(md5.as_str())))
            .collect();
        for _ in 0..20 {
            let err = open(&path, 8)
                .check(&header(&contigs), Path::new("in.cram"))
                .unwrap_err()
                .to_string();
            assert!(
                err.contains("contig `c0` has 4 bp, not its indexed 8"),
                "{err}"
            );
        }
    }

    #[test]
    fn test_each_contig_is_hashed_once_across_checks() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">a\nACGT\n>b\nTT\n");
        let check = open(&path, 2);
        let md5 = to_hex(&Md5::digest(b"ACGT").into());
        let both = header(&[("a", 4, Some(&md5)), ("b", 2, None)]);
        check.check(&both, Path::new("1.cram")).unwrap();
        assert!(check.digests[0].get().is_some());
        assert!(check.digests[1].get().is_none());
        std::fs::remove_file(&path).unwrap();
        let hashed_only = header(&[("a", 4, Some(&md5))]);
        check.check(&hashed_only, Path::new("2.cram")).unwrap();
    }

    #[test]
    fn test_parse_md5() {
        let digest = parse_md5(FIXTURE_MD5.as_bytes()).unwrap();
        assert_eq!(to_hex(&digest), FIXTURE_MD5);
        assert!(parse_md5(b"zz82f39b96c07c57c14cfb6ec9e0caa0").is_none());
        assert!(parse_md5(b"e1a8").is_none());
    }
}
