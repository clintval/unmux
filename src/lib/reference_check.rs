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

/// Checks CRAM headers against an indexed FASTA, hashing each contig at most once on `threads` threads.
pub(crate) struct HeaderCheck {
    path: PathBuf,
    index: fai::Index,
    digests: Vec<OnceLock<Md5Digest>>,
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
        let expected = expected_contigs(header, input)?;
        let positions = self.positions_for(&expected, input)?;
        let hashed: Vec<(&ExpectedContig, usize)> = expected
            .iter()
            .zip(positions)
            .filter(|(contig, _)| contig.md5.is_some())
            .collect();
        let unhashed = expected.len() - hashed.len();
        if unhashed > 0 {
            log::warn!(
                "{unhashed} of {} @SQ contig(s) in {} have no M5, so only their lengths were checked against the --reference FASTA",
                expected.len(),
                input.display()
            );
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .context("failed to build the reference check thread pool")?;
        let digests = pool.install(|| {
            hashed
                .par_iter()
                .with_max_len(1)
                .map(|&(_, i)| self.cached_digest(i))
                .collect::<Result<Vec<Md5Digest>>>()
        })?;
        for ((contig, _), actual) in hashed.iter().zip(&digests) {
            if let Some(md5) = contig.md5.filter(|md5| md5 != actual) {
                bail!(
                    "contig `{}` has MD5 {} in the --reference FASTA {} but @SQ M5:{} in {}; the FASTA is not the one the CRAM was written against",
                    contig.name,
                    to_hex(actual),
                    self.path.display(),
                    to_hex(&md5),
                    input.display()
                );
            }
        }
        Ok(())
    }

    /// The digest of the `i`th `.fai` record, hashing it only on first use.
    fn cached_digest(&self, i: usize) -> Result<Md5Digest> {
        if let Some(digest) = self.digests[i].get() {
            return Ok(*digest);
        }
        let digest = self.digest(&self.index.as_ref()[i])?;
        Ok(*self.digests[i].get_or_init(|| digest))
    }

    /// The `.fai` position of each contig, failing on one missing or of another length.
    fn positions_for(&self, expected: &[ExpectedContig], input: &Path) -> Result<Vec<usize>> {
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
                        input.display(),
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
                        input.display()
                    );
                }
                Ok(i)
            })
            .collect()
    }

    /// The MD5 of one sequence as htslib computes `M5`: uppercase, bytes outside `!`..=`~` dropped.
    fn digest(&self, record: &fai::Record) -> Result<Md5Digest> {
        let open_context = || format!("failed to open {}", self.path.display());
        let mut reader = match self.path.extension().and_then(|ext| ext.to_str()) {
            Some("gz" | "bgz") => fasta::io::indexed_reader::Builder::default()
                .set_index(fai::Index::default())
                .build_from_path(&self.path)
                .with_context(open_context)?
                .into_inner(),
            _ => fasta::io::BufReader::Uncompressed(BufReader::with_capacity(
                READ_BUFFER,
                File::open(&self.path).with_context(open_context)?,
            )),
        };
        reader
            .seek(SeekFrom::Start(record.offset()))
            .with_context(|| format!("failed to seek in {}", self.path.display()))?;
        let mut remaining = record.length();
        let mut hasher = Md5::new();
        let mut normalized = Vec::new();
        while remaining > 0 {
            let buf = reader
                .fill_buf()
                .with_context(|| format!("failed to read {}", self.path.display()))?;
            if buf.is_empty() {
                bail!(
                    "contig `{}` ends {remaining} bp early in {}; its .fai index is stale",
                    record.name(),
                    self.path.display()
                );
            }
            normalized.clear();
            let mut consumed = 0;
            while consumed < buf.len() && remaining > 0 {
                let rest = &buf[consumed..];
                let run = rest
                    .iter()
                    .position(|byte| !byte.is_ascii_graphic())
                    .unwrap_or(rest.len());
                let take = run.min(usize::try_from(remaining).unwrap_or(usize::MAX));
                normalized.extend_from_slice(&rest[..take]);
                remaining -= take as u64;
                consumed += if take == run && run < rest.len() {
                    take + 1
                } else {
                    take
                };
            }
            normalized.make_ascii_uppercase();
            hasher.update(&normalized);
            reader.consume(consumed);
        }
        Ok(hasher.finalize().into())
    }
}

/// One `@SQ` line's name, length and `M5`.
struct ExpectedContig {
    name: String,
    length: u64,
    md5: Option<Md5Digest>,
}

/// The `@SQ` lines of `header`, in header order.
fn expected_contigs(header: &sam::Header, input: &Path) -> Result<Vec<ExpectedContig>> {
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
                            input.display()
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
    fn test_header_without_sq_passes() {
        check(&[]).unwrap();
    }

    #[test]
    fn test_digest_matches_htslib_normalization() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">a\nacgT\r\nNN\r\n>b desc\nAC\n");
        let check = open(&path, 1);
        let records = check.index.as_ref();
        let digest = |i: usize| to_hex(&check.digest(&records[i]).unwrap());
        assert_eq!(digest(0), to_hex(&Md5::digest(b"ACGTNN").into()));
        assert_eq!(digest(1), to_hex(&Md5::digest(b"AC").into()));
    }

    #[test]
    fn test_digest_detects_a_stale_index() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">a\nACGT\n");
        std::fs::write(&path, ">a\nAC").unwrap();
        let check = open(&path, 1);
        let err = check
            .digest(&check.index.as_ref()[0])
            .unwrap_err()
            .to_string();
        assert!(err.contains("ends 2 bp early"), "{err}");
    }

    #[test]
    fn test_each_contig_is_hashed_once_across_checks() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_fasta(dir.path(), ">a\nACGT\n>b\nTT\n");
        let check = open(&path, 2);
        let md5 = to_hex(&Md5::digest(b"ACGT").into());
        let header = header(&[("a", 4, Some(&md5)), ("b", 2, None)]);
        check.check(&header, Path::new("1.cram")).unwrap();
        assert!(check.digests[0].get().is_some());
        assert!(check.digests[1].get().is_none());
        std::fs::remove_file(&path).unwrap();
        check.check(&header, Path::new("2.cram")).unwrap();
    }

    #[test]
    fn test_parse_md5() {
        let digest = parse_md5(FIXTURE_MD5.as_bytes()).unwrap();
        assert_eq!(to_hex(&digest), FIXTURE_MD5);
        assert!(parse_md5(b"zz82f39b96c07c57c14cfb6ec9e0caa0").is_none());
        assert!(parse_md5(b"e1a8").is_none());
    }
}
