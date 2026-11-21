//! CRAM decoding against a `--reference` FASTA, holding one reference sequence at a time.

use std::fs::File;
use std::io::{self, Seek, SeekFrom};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use anyhow::{anyhow, Context, Result};
use noodles::core::Region;
use noodles::sam::alignment::record::MappingQuality;
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
}

impl fasta::repository::Adapter for BoundedAdapter {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        let Some(reference) = &self.reference else {
            return Some(Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "decoding needs reference sequence `{}`; pass the FASTA the CRAM was written \
                     against with --reference",
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
    /// A decoder reading sequences from `reference`, when given.
    pub(crate) fn new(reference: Option<Reference>) -> Self {
        let loads_left = Arc::new(AtomicUsize::new(0));
        let repository = fasta::Repository::new(BoundedAdapter {
            reference: reference.clone(),
            loads_left: Arc::clone(&loads_left),
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
                    for record in &records {
                        out.push(owned(header, record)?);
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
        let records =
            decode(pass(Some(name.clone()))).context("failed to read CRAM slice records")?;
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
                *slot = Some(owned(header, record)?);
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

/// Convert a decoded CRAM record into an owned record.
fn owned(header: &sam::Header, record: &cram::Record<'_>) -> Result<RecordBuf> {
    let mut record = RecordBuf::try_from_alignment_record(header, record)
        .context("failed to convert CRAM record")?;
    // CRAM stores no MAPQ for an unmapped record; noodles reports 255, htslib 0.
    if record.flags().is_unmapped() && record.mapping_quality().is_none() {
        *record.mapping_quality_mut() = MappingQuality::new(0);
    }
    Ok(record)
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
        let mut decoder = CramDecoder::new(Some(reference));
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

    #[test]
    fn test_reference_reads_a_whole_sequence() {
        let reference = Reference::load(&test_data("contigs.fa")).unwrap();
        let record = reference.sequence(b"ctgB").unwrap();
        assert_eq!(record.sequence().len(), 200);
        assert!(reference.sequence(b"missing").is_err());
    }
}
