//! CRAM decoding against a `--reference` FASTA, holding one reference sequence at a time.

use std::fs::File;
use std::io::{self, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use anyhow::{anyhow, Context, Result};
use noodles::core::Region;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::{Record as _, RecordBuf};
use noodles::{cram, fasta, sam};

/// An indexed FASTA (with a `.fai`) that CRAM inputs decode against; clones share one reader.
#[derive(Clone)]
pub struct Reference {
    reader: Arc<Mutex<fasta::io::IndexedReader<fasta::io::BufReader<File>>>>,
}

impl Reference {
    /// Open an indexed FASTA (`PATH.fai` alongside).
    pub fn load(path: &Path) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| {
                format!(
                    "failed to open the --reference FASTA (it needs a .fai index): {}",
                    path.display()
                )
            })?;
        Ok(Self {
            reader: Arc::new(Mutex::new(reader)),
        })
    }

    /// Read one whole sequence by name into a buffer sized from the index.
    fn sequence(&self, name: &[u8]) -> io::Result<fasta::Record> {
        let mut reader = self
            .reader
            .lock()
            .map_err(|_| io::Error::other("the --reference reader was poisoned"))?;
        let length = reader
            .index()
            .as_ref()
            .iter()
            .find(|record| AsRef::<[u8]>::as_ref(record.name()) == name)
            .map_or(0, |record| record.length());
        let position = reader.index().query(&Region::new(name, ..))?;
        reader.get_mut().seek(SeekFrom::Start(position))?;
        let mut sequence = Vec::with_capacity(usize::try_from(length).unwrap_or(0));
        reader.read_sequence(&mut sequence)?;
        Ok(fasta::Record::new(
            fasta::record::Definition::new(name, None),
            fasta::record::Sequence::from(sequence),
        ))
    }
}

/// Raised by [`BoundedAdapter`] when one slice asks for a second sequence.
#[derive(Debug)]
struct SecondSequence;

impl std::fmt::Display for SecondSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("a CRAM slice spans more than one reference sequence")
    }
}

impl std::error::Error for SecondSequence {}

fn is_second_sequence(error: &io::Error) -> bool {
    error
        .get_ref()
        .is_some_and(|inner| inner.is::<SecondSequence>())
}

/// Loads at most `loads_left` sequences, then refuses with [`SecondSequence`].
struct BoundedAdapter {
    reference: Option<Reference>,
    loads_left: Arc<AtomicUsize>,
    path: PathBuf,
}

impl fasta::repository::Adapter for BoundedAdapter {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        let Some(reference) = &self.reference else {
            return Some(Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "CRAM input {} needs reference sequence `{}` to decode; pass the FASTA it was \
                     written against with --reference",
                    self.path.display(),
                    String::from_utf8_lossy(name)
                ),
            )));
        };
        if self.loads_left.load(Ordering::Relaxed) == 0 {
            return Some(Err(io::Error::other(SecondSequence)));
        }
        self.loads_left.fetch_sub(1, Ordering::Relaxed);
        Some(reference.sequence(name))
    }
}

/// Loads `target`, answers other names with an empty sequence, and records every name asked for.
struct PassAdapter {
    reference: Reference,
    target: Option<Vec<u8>>,
    requested: Arc<Mutex<Vec<Vec<u8>>>>,
}

impl fasta::repository::Adapter for PassAdapter {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        if let Ok(mut requested) = self.requested.lock() {
            requested.push(name.to_vec());
        }
        if self.target.as_deref() == Some(name) {
            return Some(self.reference.sequence(name));
        }
        Some(Ok(fasta::Record::new(
            fasta::record::Definition::new(name, None),
            fasta::record::Sequence::default(),
        )))
    }
}

/// Decodes one CRAM input's containers into owned records.
pub(crate) struct CramDecoder {
    reference: Option<Reference>,
    repository: fasta::Repository,
    loads_left: Arc<AtomicUsize>,
}

impl CramDecoder {
    /// A decoder for the CRAM at `path`, reading sequences from `reference` when given.
    pub(crate) fn new(reference: Option<Reference>, path: &Path) -> Self {
        let loads_left = Arc::new(AtomicUsize::new(0));
        let repository = fasta::Repository::new(BoundedAdapter {
            reference: reference.clone(),
            loads_left: Arc::clone(&loads_left),
            path: path.to_path_buf(),
        });
        Self {
            reference,
            repository,
            loads_left,
        }
    }

    /// Decode every slice of `container` into owned records, in file order.
    pub(crate) fn decode(
        &mut self,
        container: &cram::io::reader::Container,
        header: &sam::Header,
    ) -> Result<Vec<RecordBuf>> {
        let compression_header = container
            .compression_header()
            .context("failed to read CRAM compression header")?;
        let mut out = Vec::new();
        for slice in container.slices() {
            let slice = slice.context("failed to read CRAM slice")?;
            let (core, external) = slice
                .decode_blocks()
                .context("failed to decode CRAM slice blocks")?;
            let decode = |repository: fasta::Repository| {
                slice.records(repository, header, &compression_header, &core, &external)
            };
            let may_load = self.repository.is_empty();
            let decoded = match self.try_decode(&decode, may_load) {
                Err(error) if is_second_sequence(&error) && !may_load => {
                    self.repository.clear();
                    self.try_decode(&decode, true)
                }
                decoded => decoded,
            };
            match decoded {
                Ok(records) => {
                    let cached = |name: &[u8]| self.repository.get(name).and_then(Result::ok);
                    for record in &records {
                        out.push(owned(header, record, &cached)?);
                    }
                }
                Err(error) if is_second_sequence(&error) => {
                    self.repository.clear();
                    let reference = self
                        .reference
                        .clone()
                        .ok_or_else(|| anyhow!("a multi-reference CRAM slice needs --reference"))?;
                    decode_per_sequence(&reference, header, decode, &mut out)?;
                }
                Err(error) => return Err(error).context("failed to read CRAM slice records"),
            }
        }
        Ok(out)
    }

    /// Decode one slice, loading a sequence only when `may_load` is set.
    fn try_decode<'c, F>(&self, decode: &F, may_load: bool) -> io::Result<Vec<cram::Record<'c>>>
    where
        F: Fn(fasta::Repository) -> io::Result<Vec<cram::Record<'c>>>,
    {
        self.loads_left
            .store(usize::from(may_load), Ordering::Relaxed);
        let decoded = decode(self.repository.clone());
        self.loads_left.store(0, Ordering::Relaxed);
        decoded
    }
}

/// Decode a multi-reference slice in one pass per sequence, so one sequence is resident at a time.
fn decode_per_sequence<'c, F>(
    reference: &Reference,
    header: &sam::Header,
    decode: F,
    out: &mut Vec<RecordBuf>,
) -> Result<()>
where
    F: Fn(fasta::Repository) -> io::Result<Vec<cram::Record<'c>>>,
{
    let requested = Arc::new(Mutex::new(Vec::new()));
    let pass = |target: Option<Vec<u8>>| {
        fasta::Repository::new(PassAdapter {
            reference: reference.clone(),
            target,
            requested: Arc::clone(&requested),
        })
    };
    let count = decode(pass(None))
        .context("failed to read CRAM slice records")?
        .len();
    let mut names: Vec<Vec<u8>> = requested
        .lock()
        .map_err(|_| anyhow!("the CRAM reference bookkeeping was poisoned"))?
        .clone();
    let mut seen = std::collections::HashSet::new();
    names.retain(|name| seen.insert(name.clone()));

    let mut slots: Vec<Option<RecordBuf>> = vec![None; count];
    for (i, name) in names.iter().enumerate() {
        let repository = pass(Some(name.clone()));
        let records = decode(repository.clone()).context("failed to read CRAM slice records")?;
        let loaded = |on: &[u8]| {
            (on == name.as_slice())
                .then(|| repository.get(on).and_then(Result::ok))
                .flatten()
        };
        for (slot, record) in slots.iter_mut().zip(&records) {
            if slot.is_some() {
                continue;
            }
            let on = sequence_name(header, record)?;
            let wanted = match on {
                Some(on) if names.iter().any(|n| n.as_slice() == on) => on == name.as_slice(),
                _ => i == 0,
            };
            if wanted {
                *slot = Some(owned(header, record, &loaded)?);
            }
        }
    }
    for slot in slots {
        out.push(slot.ok_or_else(|| anyhow!("a CRAM slice record was not decoded"))?);
    }
    Ok(())
}

/// The reference sequence name a mapped record is aligned to.
fn sequence_name<'h>(
    header: &'h sam::Header,
    record: &cram::Record<'_>,
) -> Result<Option<&'h [u8]>> {
    if record.flags()?.is_unmapped() {
        return Ok(None);
    }
    let Some(id) = record.reference_sequence_id(header).transpose()? else {
        return Ok(None);
    };
    Ok(header
        .reference_sequences()
        .get_index(id)
        .map(|(name, _)| name.as_ref()))
}

/// The loaded reference sequence a record decoded against, looked up by name.
type Loaded<'a> = &'a dyn Fn(&[u8]) -> Option<Arc<fasta::record::Sequence>>;

/// Convert a decoded CRAM record into an owned record, as htslib would decode it.
fn owned(header: &sam::Header, record: &cram::Record<'_>, loaded: Loaded<'_>) -> Result<RecordBuf> {
    let mut record = RecordBuf::try_from_alignment_record(header, record)
        .context("failed to convert CRAM record")?;
    // CRAM stores no MAPQ for an unmapped record; noodles reports 255, htslib 0.
    if record.flags().is_unmapped() && record.mapping_quality().is_none() {
        *record.mapping_quality_mut() = MappingQuality::new(0);
    }
    record.sequence_mut().as_mut().make_ascii_uppercase();
    regenerate_md_nm(header, &mut record, loaded);
    Ok(record)
}

/// Add the MD and NM tags htslib regenerates when it decodes a mapped record against its reference.
fn regenerate_md_nm(header: &sam::Header, record: &mut RecordBuf, loaded: Loaded<'_>) {
    let skip = take_cram_flags(record);
    let want_md = skip & 1 == 0 && record.data().get(&Tag::MISMATCHED_POSITIONS).is_none();
    let want_nm = skip & 2 == 0 && record.data().get(&Tag::EDIT_DISTANCE).is_none();
    if !(want_md || want_nm) || record.flags().is_unmapped() || record.sequence().is_empty() {
        return;
    }
    let Some(sequence) = record
        .reference_sequence_id()
        .and_then(|id| header.reference_sequences().get_index(id))
        .and_then(|(name, _)| loaded(name))
    else {
        return;
    };
    if let Some((md, nm)) = md_nm(record, (*sequence).as_ref()) {
        add_md_nm(header, record, want_md.then_some(md), want_nm.then_some(nm));
    }
}

/// Remove htslib's private `cF` tag, returning its bits: 1 keeps MD absent, 2 keeps NM absent.
fn take_cram_flags(record: &mut RecordBuf) -> i64 {
    record
        .data_mut()
        .remove(&Tag::new(b'c', b'F'))
        .and_then(|(_, value)| value.as_int())
        .unwrap_or(0)
}

/// The MD string and edit distance htslib derives from a record's CIGAR, SEQ and reference.
fn md_nm(record: &RecordBuf, reference: &[u8]) -> Option<(Vec<u8>, u32)> {
    let read = record.sequence().as_ref();
    let mut ref_pos = usize::from(record.alignment_start()?) - 1;
    let mut read_pos = 0;
    let mut md = Vec::new();
    let mut nm = 0u32;
    let mut run = 0u32;
    for op in record.cigar().as_ref() {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..len {
                    let base = reference.get(ref_pos)?.to_ascii_uppercase();
                    if read.get(read_pos)?.to_ascii_uppercase() != base || base == b'N' {
                        md.extend_from_slice(run.to_string().as_bytes());
                        md.push(base);
                        run = 0;
                        nm += 1;
                    } else {
                        run += 1;
                    }
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Kind::Deletion => {
                let bases = reference.get(ref_pos..ref_pos + len)?;
                md.extend_from_slice(run.to_string().as_bytes());
                md.push(b'^');
                md.extend(bases.iter().map(u8::to_ascii_uppercase));
                run = 0;
                nm += u32::try_from(len).ok()?;
                ref_pos += len;
            }
            Kind::Insertion => {
                nm += u32::try_from(len).ok()?;
                read_pos += len;
            }
            Kind::SoftClip => read_pos += len,
            Kind::Skip => ref_pos += len,
            Kind::HardClip | Kind::Pad => {}
        }
    }
    md.extend_from_slice(run.to_string().as_bytes());
    Some((md, nm))
}

/// Add MD and NM where htslib does: after the stored tags, ahead of a declared `RG:Z`.
fn add_md_nm(header: &sam::Header, record: &mut RecordBuf, md: Option<Vec<u8>>, nm: Option<u32>) {
    let mut fields: Vec<(Tag, Value)> = record
        .data()
        .iter()
        .map(|(tag, value)| (tag, value.clone()))
        .collect();
    let at = match fields.last() {
        Some((tag, Value::String(id)))
            if *tag == Tag::READ_GROUP && header.read_groups().contains_key(id.as_slice()) =>
        {
            fields.len() - 1
        }
        _ => fields.len(),
    };
    let md = md.map(|md| (Tag::MISMATCHED_POSITIONS, Value::String(md.into())));
    let nm = nm.map(|nm| {
        let value = match (u8::try_from(nm), u16::try_from(nm)) {
            (Ok(nm), _) => Value::UInt8(nm),
            (_, Ok(nm)) => Value::UInt16(nm),
            _ => Value::UInt32(nm),
        };
        (Tag::EDIT_DISTANCE, value)
    });
    fields.splice(at..at, md.into_iter().chain(nm));
    *record.data_mut() = fields.into_iter().collect();
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn test_data(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name)
    }

    /// Decode a test CRAM, returning its records and the most sequences held after any container.
    fn decode_all(name: &str) -> (Vec<RecordBuf>, usize) {
        let reference = Reference::load(&test_data("contigs.fa")).unwrap();
        let mut reader = cram::io::Reader::new(File::open(test_data(name)).unwrap());
        let header = reader.read_header().unwrap();
        let mut decoder = CramDecoder::new(Some(reference), &test_data(name));
        let mut container = cram::io::reader::Container::default();
        let mut records = Vec::new();
        let mut most = 0;
        while reader.read_container(&mut container).unwrap() != 0 {
            records.extend(decoder.decode(&container, &header).unwrap());
            most = most.max(decoder.repository.len());
        }
        (records, most)
    }

    #[test]
    fn test_single_reference_slices_keep_one_sequence() {
        let (records, most) = decode_all("contigs.cram");
        assert_eq!(records.len(), 29);
        assert_eq!(most, 1);
    }

    #[test]
    fn test_multi_reference_slice_decodes_like_single_reference_slices() {
        let (single, _) = decode_all("contigs.cram");
        let (multi, most) = decode_all("contigs.multi-ref.cram");
        assert_eq!(multi, single);
        assert_eq!(most, 0);
    }

    /// The reference behind the MD/NM tests: an N run at 5-6 and lowercase bases at 11-14.
    const CONTIG: &[u8] = b"ACGTNNACGTacgtACGTACGTACGTACGT";

    /// Parse one SAM record against a header declaring contig `c` and read group `rg1`.
    fn parse(line: &str) -> (sam::Header, RecordBuf) {
        let text = format!("@SQ\tSN:c\tLN:30\n@RG\tID:rg1\n{line}\n");
        let mut reader = sam::io::Reader::new(text.as_bytes());
        let header = reader.read_header().unwrap();
        let mut record = RecordBuf::default();
        reader.read_record_buf(&header, &mut record).unwrap();
        (header, record)
    }

    /// Regenerate MD/NM on a SAM record and render its tags as SAM text.
    fn regenerated(line: &str) -> String {
        use noodles::sam::alignment::io::Write as _;
        let (header, mut record) = parse(line);
        let contig = Arc::new(fasta::record::Sequence::from(CONTIG.to_vec()));
        regenerate_md_nm(&header, &mut record, &|_: &[u8]| Some(Arc::clone(&contig)));
        let mut out = Vec::new();
        sam::io::Writer::new(&mut out)
            .write_alignment_record(&header, &record)
            .unwrap();
        let text = String::from_utf8(out).unwrap();
        text.trim_end()
            .split('\t')
            .skip(11)
            .collect::<Vec<_>>()
            .join(" ")
    }

    #[test]
    fn test_md_nm_count_reference_n_as_a_mismatch() {
        assert_eq!(
            regenerated("r\t0\tc\t3\t60\t6M\t*\t0\t0\tGTNNAC\t*\tRG:Z:rg1"),
            "MD:Z:2N0N2 NM:i:2 RG:Z:rg1"
        );
        assert_eq!(
            regenerated("r\t0\tc\t1\t60\t6M\t*\t0\t0\tANGTAA\t*"),
            "MD:Z:1C2N0N0 NM:i:3"
        );
    }

    #[test]
    fn test_md_nm_follow_indels_clips_and_skips() {
        assert_eq!(
            regenerated("r\t0\tc\t1\t60\t2S3M2I2M2D2M3N2M\t*\t0\t0\tTTACGGGTNCGGT\t*"),
            "MD:Z:4N0^NA4 NM:i:5"
        );
        assert_eq!(
            regenerated("r\t0\tc\t7\t60\t3H4M\t*\t0\t0\tACGT\t*"),
            "MD:Z:4 NM:i:0"
        );
    }

    #[test]
    fn test_md_nm_compare_against_an_uppercased_reference() {
        assert_eq!(
            regenerated("r\t0\tc\t9\t60\t6M\t*\t0\t0\tGTACTA\t*"),
            "MD:Z:4G0T0 NM:i:2"
        );
        assert_eq!(
            regenerated("r\t0\tc\t10\t60\t3D2M\t*\t0\t0\tGT\t*"),
            "MD:Z:0^TAC2 NM:i:3"
        );
    }

    #[test]
    fn test_md_nm_skip_unmapped_and_sequenceless_records() {
        assert_eq!(regenerated("r\t4\tc\t1\t0\t*\t*\t0\t0\tACGT\t*"), "");
        assert_eq!(regenerated("r\t256\tc\t1\t0\t4M\t*\t0\t0\t*\t*"), "");
    }

    #[test]
    fn test_md_nm_keep_stored_tags_and_add_the_missing_one() {
        assert_eq!(
            regenerated("r\t0\tc\t1\t60\t4M\t*\t0\t0\tACGA\t*\tMD:Z:9\tRG:Z:rg1"),
            "MD:Z:9 NM:i:1 RG:Z:rg1"
        );
        assert_eq!(
            regenerated("r\t0\tc\t1\t60\t4M\t*\t0\t0\tACGA\t*\tNM:i:7"),
            "NM:i:7 MD:Z:3T0"
        );
    }

    #[test]
    fn test_md_nm_follow_an_undeclared_read_group() {
        assert_eq!(
            regenerated("r\t0\tc\t1\t60\t4M\t*\t0\t0\tACGT\t*\tXA:Z:x\tRG:Z:ghost"),
            "XA:Z:x RG:Z:ghost MD:Z:4 NM:i:0"
        );
    }

    #[test]
    fn test_cram_flags_tag_suppresses_md_or_nm_and_is_removed() {
        let line = "r\t0\tc\t1\t60\t4M\t*\t0\t0\tACGT\t*\tcF:i:1\tRG:Z:rg1";
        assert_eq!(regenerated(line), "NM:i:0 RG:Z:rg1");
        let line = "r\t0\tc\t1\t60\t4M\t*\t0\t0\tACGT\t*\tcF:i:3";
        assert_eq!(regenerated(line), "");
    }

    #[test]
    fn test_nm_takes_the_narrowest_unsigned_type() {
        let (header, mut record) = parse("r\t0\tc\t1\t60\t4M\t*\t0\t0\tACGT\t*");
        for (nm, want) in [
            (255, Value::UInt8(255)),
            (256, Value::UInt16(256)),
            (70_000, Value::UInt32(70_000)),
        ] {
            record.data_mut().clear();
            add_md_nm(&header, &mut record, None, Some(nm));
            assert_eq!(record.data().get(&Tag::EDIT_DISTANCE), Some(&want));
        }
    }

    #[test]
    fn test_reference_reads_a_whole_sequence() {
        let reference = Reference::load(&test_data("contigs.fa")).unwrap();
        let record = reference.sequence(b"ctgB").unwrap();
        assert_eq!(record.sequence().len(), 200);
        assert!(reference.sequence(b"missing").is_err());
    }
}
