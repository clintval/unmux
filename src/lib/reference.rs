//! CRAM decoding against a `--reference` FASTA, holding only the reference sequences slices need.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, Seek, SeekFrom};
use std::path::{Path, PathBuf};
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
    fn sequence(&self, name: &[u8]) -> io::Result<fasta::record::Sequence> {
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
        Ok(fasta::record::Sequence::from(sequence))
    }
}

/// The reference sequences a [`CramDecoder`] holds, by name, with their lengths.
type Held = Arc<Mutex<HashMap<Vec<u8>, usize>>>;

/// Loads sequences from the `--reference` FASTA, or takes back ones carried over a rebuild.
struct LoadingAdapter {
    reference: Option<Reference>,
    path: PathBuf,
    lengths: Arc<HashMap<Vec<u8>, usize>>,
    carried: HashMap<Vec<u8>, Arc<fasta::record::Sequence>>,
    held: Held,
}

impl fasta::repository::Adapter for LoadingAdapter {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        let sequence = match (self.carried.remove(name), &self.reference) {
            (Some(sequence), _) => Arc::unwrap_or_clone(sequence),
            (None, Some(reference)) => match reference.sequence(name) {
                Ok(sequence) => match self.lengths.get(name) {
                    Some(&length) if sequence.len() < length => {
                        return Some(Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!(
                                "reference sequence `{}` is {} bases in the --reference FASTA but \
                                 {} in the @SQ header of CRAM input {}; pass the FASTA it was \
                                 written against",
                                String::from_utf8_lossy(name),
                                sequence.len(),
                                length,
                                self.path.display()
                            ),
                        )))
                    }
                    _ => sequence,
                },
                Err(error) => return Some(Err(error)),
            },
            (None, None) => {
                return Some(Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "CRAM input {} needs reference sequence `{}` to decode; pass the FASTA it \
                         was written against with --reference",
                        self.path.display(),
                        String::from_utf8_lossy(name)
                    ),
                )))
            }
        };
        if let Ok(mut held) = self.held.lock() {
            held.insert(name.to_vec(), sequence.len());
        }
        Some(Ok(fasta::Record::new(
            fasta::record::Definition::new(name, None),
            sequence,
        )))
    }
}

/// Serves padded sequences, and an empty stand-in for any other name.
struct PaddedAdapter(HashMap<Vec<u8>, fasta::record::Sequence>);

impl fasta::repository::Adapter for PaddedAdapter {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        Some(Ok(fasta::Record::new(
            fasta::record::Definition::new(name, None),
            self.0.remove(name).unwrap_or_default(),
        )))
    }
}

/// Serves zeroed stand-ins of the given lengths for slices that need no reference: noodles still
/// asks for the sequence of every record in a multi-reference slice and indexes it by position.
struct BlankAdapter(HashMap<Vec<u8>, usize>);

impl fasta::repository::Adapter for BlankAdapter {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        let length = self.0.get(name).copied().unwrap_or(0);
        Some(Ok(fasta::Record::new(
            fasta::record::Definition::new(name, None),
            fasta::record::Sequence::from(vec![0; length]),
        )))
    }
}

/// Raised by [`ProbeAdapter`] with the first sequence a slice asks for.
#[derive(Debug)]
struct Requested(Vec<u8>);

impl std::fmt::Display for Requested {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("a CRAM slice asked for a reference sequence")
    }
}

impl std::error::Error for Requested {}

/// Refuses every sequence, so a slice decodes only if it needs no reference.
struct ProbeAdapter;

impl fasta::repository::Adapter for ProbeAdapter {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        Some(Err(io::Error::other(Requested(name.to_vec()))))
    }
}

/// The sequence a [`ProbeAdapter`] refused, if that is what failed a decode.
fn requested(error: &io::Error) -> Option<&[u8]> {
    error
        .get_ref()
        .and_then(|inner| inner.downcast_ref::<Requested>())
        .map(|requested| requested.0.as_slice())
}

/// Decodes one CRAM input's containers into owned records.
pub(crate) struct CramDecoder {
    reference: Option<Reference>,
    path: PathBuf,
    lengths: Arc<HashMap<Vec<u8>, usize>>,
    repository: fasta::Repository,
    held: Held,
    blank: Option<fasta::Repository>,
    blank_lengths: HashMap<Vec<u8>, usize>,
}

impl CramDecoder {
    /// A decoder for the CRAM at `path` with `header`, reading sequences from `reference` when
    /// given.
    pub(crate) fn new(reference: Option<Reference>, path: &Path, header: &sam::Header) -> Self {
        let lengths: HashMap<Vec<u8>, usize> = header
            .reference_sequences()
            .iter()
            .map(|(name, sequence)| (name.to_vec(), sequence.length().get()))
            .collect();
        let mut decoder = Self {
            reference,
            path: path.to_path_buf(),
            lengths: Arc::new(lengths.clone()),
            repository: fasta::Repository::default(),
            held: Held::default(),
            blank: None,
            blank_lengths: lengths,
        };
        decoder.repository = decoder.loading_repository(HashMap::new());
        decoder
    }

    fn loading_repository(
        &self,
        carried: HashMap<Vec<u8>, Arc<fasta::record::Sequence>>,
    ) -> fasta::Repository {
        fasta::Repository::new(LoadingAdapter {
            reference: self.reference.clone(),
            path: self.path.clone(),
            lengths: Arc::clone(&self.lengths),
            carried,
            held: Arc::clone(&self.held),
        })
    }

    fn held(&self) -> Result<std::sync::MutexGuard<'_, HashMap<Vec<u8>, usize>>> {
        self.held
            .lock()
            .map_err(|_| anyhow!("the CRAM reference bookkeeping was poisoned"))
    }

    /// A held sequence, never loading one.
    fn held_sequence(&self, name: &[u8]) -> Option<Arc<fasta::record::Sequence>> {
        if !self.held().ok()?.contains_key(name) {
            return None;
        }
        self.repository.get(name).and_then(Result::ok)
    }

    /// Decode every slice of `container` into owned records, in file order.
    ///
    /// Each slice decodes once against every sequence it references. A probe with no sequences
    /// first finds the first one it needs, so a single-reference container, or sorted input, drops
    /// what the slice cannot use before loading more; afterwards the unused sequences are dropped.
    /// A container written without the reference never loads one.
    pub(crate) fn decode(
        &mut self,
        container: &cram::io::reader::Container,
        header: &sam::Header,
    ) -> Result<Vec<RecordBuf>> {
        let compression_header = container
            .compression_header()
            .context("failed to read CRAM compression header")?;
        let regenerate = uses_external_reference(&compression_header)?;
        let single_reference = !container.header().reference_sequence_context().is_many();
        let sorted = is_coordinate_sorted(header);
        let mut out = Vec::new();
        for slice in container.slices() {
            let slice = slice.context("failed to read CRAM slice")?;
            let (core, external) = slice
                .decode_blocks()
                .context("failed to decode CRAM slice blocks")?;
            let decode = |repository: fasta::Repository| {
                slice.records(repository, header, &compression_header, &core, &external)
            };
            let records = if !regenerate {
                let records = decode(self.blank_repository(&[]))
                    .context("failed to read CRAM slice records")?;
                let ends: Vec<(&[u8], usize)> = records
                    .iter()
                    .map(|record| aligned_end(header, record))
                    .filter_map(Result::transpose)
                    .collect::<Result<_>>()?;
                if sorted {
                    self.blank = None;
                }
                if ends
                    .iter()
                    .all(|(name, end)| self.blank_lengths.get(*name) >= Some(end))
                {
                    Ok(records)
                } else {
                    drop(records);
                    decode(self.blank_repository(&ends))
                }
            } else {
                match decode(fasta::Repository::new(ProbeAdapter)) {
                    Err(error) => match requested(&error) {
                        Some(name) => {
                            if single_reference {
                                self.retain(&HashSet::from([name]))?;
                            } else if sorted {
                                self.retain_from(header, name)?;
                            }
                            decode(self.repository.clone())
                        }
                        None => Err(error),
                    },
                    records => records,
                }
            }
            .context("failed to read CRAM slice records")?;
            let lengths = self.held()?.clone();
            let mut padded: HashMap<&[u8], usize> = HashMap::new();
            for record in &records {
                if let Some((name, end)) = aligned_end(header, record)? {
                    if lengths.get(name).is_some_and(|&length| end > length) {
                        let length = padded.entry(name).or_default();
                        *length = (*length).max(end);
                    }
                }
            }
            let repadded = if padded.is_empty() {
                Vec::new()
            } else {
                decode(self.padded_repository(&padded)?)
                    .context("failed to read CRAM slice records")?
            };
            let mut used = HashSet::new();
            for (i, record) in records.iter().enumerate() {
                let name = sequence_name(header, record)?;
                used.extend(name);
                let record = match name.filter(|name| padded.contains_key(name)) {
                    Some(_) => &repadded[i],
                    None => record,
                };
                let loaded = |name: &[u8]| self.held_sequence(name).filter(|_| regenerate);
                let converted = match aligned_end(header, record)? {
                    Some((name, end))
                        if end > header_length(header, name) && !padded.contains_key(name) =>
                    {
                        convert_guarded(header, record, &loaded, || past_end(record, name))
                    }
                    _ if !regenerate => {
                        convert_guarded(header, record, &loaded, || self.unreferenced(record))
                    }
                    _ => owned(header, record, &loaded),
                }?;
                // A zeroed stand-in base means the record needed reference bases its slice lacks.
                if !regenerate && converted.sequence().as_ref().contains(&0) {
                    return Err(self.unreferenced(record));
                }
                out.push(converted);
            }
            drop((records, repadded));
            self.retain(&used)?;
        }
        Ok(out)
    }

    /// The error for a record that needs reference bases its slice does not carry.
    fn unreferenced(&self, record: &cram::Record<'_>) -> anyhow::Error {
        anyhow!(
            "CRAM record `{}` in {} needs reference bases its slice does not carry",
            record_name(record),
            self.path.display()
        )
    }

    /// The zeroed stand-ins for slices that need no reference, each as long as its header length
    /// or the farthest of `ends` on it.
    fn blank_repository(&mut self, ends: &[(&[u8], usize)]) -> fasta::Repository {
        for &(name, end) in ends {
            let length = self.blank_lengths.entry(name.to_vec()).or_default();
            if end > *length {
                *length = end;
                self.blank = None;
            }
        }
        self.blank
            .get_or_insert_with(|| fasta::Repository::new(BlankAdapter(self.blank_lengths.clone())))
            .clone()
    }

    /// Held sequences extended with N to the given lengths, and an empty stand-in for any other.
    fn padded_repository(&self, padded: &HashMap<&[u8], usize>) -> Result<fasta::Repository> {
        let mut sequences = HashMap::new();
        for (&name, &length) in padded {
            let held = self
                .held_sequence(name)
                .ok_or_else(|| anyhow!("a padded CRAM reference sequence was not held"))?;
            let mut bases = Vec::with_capacity(length);
            bases.extend_from_slice((*held).as_ref());
            bases.resize(length.max(bases.len()), b'N');
            sequences.insert(name.to_vec(), fasta::record::Sequence::from(bases));
        }
        Ok(fasta::Repository::new(PaddedAdapter(sequences)))
    }

    /// Drop the held sequences ahead of `first` in header order, which sorted input never revisits.
    fn retain_from(&mut self, header: &sam::Header, first: &[u8]) -> Result<()> {
        let sequences = header.reference_sequences();
        let first = sequences.get_index_of(first);
        let needed: Vec<Vec<u8>> = self
            .held()?
            .keys()
            .filter(|name| sequences.get_index_of(name.as_slice()) >= first)
            .cloned()
            .collect();
        self.retain(&needed.iter().map(Vec::as_slice).collect())
    }

    /// Drop the held sequences outside `needed`, unless they free less than keeping the rest copies.
    fn retain(&mut self, needed: &HashSet<&[u8]>) -> Result<()> {
        let held = self.held()?.clone();
        if held.keys().all(|name| needed.contains(name.as_slice())) {
            return Ok(());
        }
        let (kept, dropped): (Vec<_>, Vec<_>) = held
            .iter()
            .partition(|(name, _)| needed.contains(name.as_slice()));
        let bases = |part: &[(&Vec<u8>, &usize)]| part.iter().map(|(_, len)| **len).sum::<usize>();
        if bases(&dropped) < bases(&kept) {
            return Ok(());
        }
        let carried: HashMap<Vec<u8>, Arc<fasta::record::Sequence>> = kept
            .iter()
            .filter_map(|(name, _)| Some(((*name).clone(), self.repository.get(name)?.ok()?)))
            .collect();
        let names: Vec<Vec<u8>> = carried.keys().cloned().collect();
        self.held()?.clear();
        self.repository = fasta::Repository::default();
        self.repository = self.loading_repository(carried);
        for name in names {
            self.repository
                .get(&name)
                .transpose()
                .context("failed to keep a CRAM reference sequence")?;
        }
        Ok(())
    }
}

/// Whether a container's slices decode against the external reference: the compression header's
/// `RR` flag, which noodles parses but keeps private, so it is read from the header's debug form.
fn uses_external_reference(
    compression_header: &cram::container::CompressionHeader,
) -> Result<bool> {
    const FLAG: &str = "external_reference_sequence_is_required: ";
    let debug = format!("{compression_header:?}");
    let at = debug
        .find(FLAG)
        .ok_or_else(|| anyhow!("failed to read the reference flag of a CRAM compression header"))?;
    Ok(debug[at + FLAG.len()..].starts_with("true"))
}

/// Whether the header declares the records sorted by coordinate.
fn is_coordinate_sorted(header: &sam::Header) -> bool {
    use noodles::sam::header::record::value::map::header::{sort_order, tag};
    header
        .header()
        .and_then(|map| map.other_fields().get(&tag::SORT_ORDER))
        .is_some_and(|order| order.as_slice() == sort_order::COORDINATE)
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

/// The sequence a mapped record is aligned to and its alignment end.
fn aligned_end<'h>(
    header: &'h sam::Header,
    record: &cram::Record<'_>,
) -> Result<Option<(&'h [u8], usize)>> {
    let Some(name) = sequence_name(header, record)? else {
        return Ok(None);
    };
    let end = record.alignment_end().transpose()?;
    Ok(end.map(|end| (name, usize::from(end))))
}

/// The length the header gives sequence `name`.
fn header_length(header: &sam::Header, name: &[u8]) -> usize {
    header
        .reference_sequences()
        .get(name)
        .map_or(usize::MAX, |sequence| sequence.length().get())
}

/// Convert a record whose decode can panic inside noodles, failing with `error` instead.
fn convert_guarded(
    header: &sam::Header,
    record: &cram::Record<'_>,
    loaded: Loaded<'_>,
    error: impl FnOnce() -> anyhow::Error,
) -> Result<RecordBuf> {
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        owned(header, record, loaded)
    }))
    .unwrap_or_else(|_| Err(error()))
}

/// The error for a record that runs past the end of `name` in a slice that embeds its reference.
fn past_end(record: &cram::Record<'_>, name: &[u8]) -> anyhow::Error {
    anyhow!(
        "CRAM record `{}` runs past the end of reference sequence `{}` in a slice that embeds its \
         reference, which unmux cannot decode; convert the input to BAM first (e.g. `samtools \
         view -b`)",
        record_name(record),
        String::from_utf8_lossy(name)
    )
}

/// A record's name, for messages.
fn record_name(record: &cram::Record<'_>) -> String {
    record.name().map(ToString::to_string).unwrap_or_default()
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

/// Add the MD and NM tags htslib regenerates when it decodes a mapped record against the external
/// reference, which `loaded` holds.
fn regenerate_md_nm(header: &sam::Header, record: &mut RecordBuf, loaded: Loaded<'_>) {
    let skip = take_cram_flags(record);
    let want_md = skip & 1 == 0 && record.data().get(&Tag::MISMATCHED_POSITIONS).is_none();
    let want_nm = skip & 2 == 0 && record.data().get(&Tag::EDIT_DISTANCE).is_none();
    if !(want_md || want_nm) || record.flags().is_unmapped() || record.sequence().is_empty() {
        return;
    }
    let Some((sequence, length)) = record
        .reference_sequence_id()
        .and_then(|id| header.reference_sequences().get_index(id))
        .and_then(|(name, map)| Some((loaded(name)?, map.length().get())))
    else {
        return;
    };
    let sequence: &[u8] = (*sequence).as_ref();
    if let Some((md, nm)) = md_nm(record, &sequence[..length.min(sequence.len())]) {
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

/// The MD string and edit distance htslib derives from a record's CIGAR, SEQ and reference; like
/// htslib, the MD string stops at the first base past the end of the reference.
fn md_nm(record: &RecordBuf, reference: &[u8]) -> Option<(Vec<u8>, u32)> {
    let read = record.sequence().as_ref();
    let mut ref_pos = usize::from(record.alignment_start()?) - 1;
    let mut read_pos = 0;
    let mut md = Vec::new();
    let mut nm = 0u32;
    let mut run = 0u32;
    let mut past_end = false;
    let end_run = |md: &mut Vec<u8>, run: &mut u32| {
        md.extend_from_slice(run.to_string().as_bytes());
        *run = 0;
    };
    for op in record.cigar().as_ref() {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..len {
                    let read_base = read.get(read_pos)?.to_ascii_uppercase();
                    match reference.get(ref_pos).map(u8::to_ascii_uppercase) {
                        _ if past_end => {}
                        None => {
                            end_run(&mut md, &mut run);
                            past_end = true;
                        }
                        Some(base) if read_base != base || base == b'N' => {
                            end_run(&mut md, &mut run);
                            md.push(base);
                            nm += 1;
                        }
                        Some(_) => run += 1,
                    }
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Kind::Deletion => {
                if !past_end {
                    end_run(&mut md, &mut run);
                    if let Some(bases) = reference.get(ref_pos..reference.len().min(ref_pos + len))
                    {
                        md.push(b'^');
                        md.extend(bases.iter().map(u8::to_ascii_uppercase));
                        nm += u32::try_from(bases.len()).ok()?;
                        past_end = bases.len() < len;
                        if past_end {
                            md.push(b'0');
                        }
                    } else {
                        past_end = true;
                    }
                }
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
    if !past_end {
        end_run(&mut md, &mut run);
    }
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
        let mut decoder = CramDecoder::new(Some(reference), &test_data(name), &header);
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
    fn test_the_reference_flag_is_read_from_the_compression_header() {
        for (name, external) in [
            ("contigs.cram", true),
            ("contigs.multi-ref.cram", true),
            ("contigs.no-ref.cram", false),
            ("contigs.embedded-ref.cram", false),
        ] {
            let mut reader = cram::io::Reader::new(File::open(test_data(name)).unwrap());
            reader.read_header().unwrap();
            let mut container = cram::io::reader::Container::default();
            reader.read_container(&mut container).unwrap();
            let compression_header = container.compression_header().unwrap();
            assert_eq!(
                uses_external_reference(&compression_header).unwrap(),
                external,
                "{name}"
            );
        }
    }

    #[test]
    fn test_slices_without_the_external_reference_get_no_md_or_nm() {
        let (records, _) = decode_all("contigs.mixed-ref.cram");
        let with_md = |records: &[RecordBuf]| {
            records
                .iter()
                .filter(|record| record.data().get(&Tag::MISMATCHED_POSITIONS).is_some())
                .count()
        };
        let parts: Vec<&[RecordBuf]> = records.chunks(29).collect();
        let md: Vec<usize> = parts.iter().map(|part| with_md(part)).collect();
        assert_eq!(md, [26, 1, 26, 1]);
    }

    #[test]
    fn test_multi_reference_slice_written_without_a_reference_never_loads_one() {
        let path = test_data("contigs.no-ref.multi-ref.cram");
        for reference in [
            None,
            Some(Reference::load(&test_data("contigs.fa")).unwrap()),
        ] {
            let mut reader = cram::io::Reader::new(File::open(&path).unwrap());
            let header = reader.read_header().unwrap();
            let mut decoder = CramDecoder::new(reference, &path, &header);
            let mut container = cram::io::reader::Container::default();
            let mut records = Vec::new();
            while reader.read_container(&mut container).unwrap() != 0 {
                records.extend(decoder.decode(&container, &header).unwrap());
            }
            assert_eq!(records.len(), 29);
            assert!(decoder.repository.is_empty());
            assert!(records
                .iter()
                .all(|record| !record.sequence().as_ref().contains(&0)));
        }
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
        assert_eq!(most, 3);
    }

    #[test]
    fn test_retain_drops_unneeded_sequences_unless_keeping_the_rest_costs_more() {
        let reference = Reference::load(&test_data("contigs.fa")).unwrap();
        let path = test_data("contigs.cram");
        let header = cram::io::Reader::new(File::open(&path).unwrap())
            .read_header()
            .unwrap();
        let mut decoder = CramDecoder::new(Some(reference), &path, &header);
        let load = |decoder: &CramDecoder, name: &[u8]| {
            decoder.repository.get(name).unwrap().unwrap();
        };
        let held = |decoder: &CramDecoder| {
            let mut names: Vec<Vec<u8>> = decoder.held().unwrap().keys().cloned().collect();
            names.sort();
            names
        };
        load(&decoder, b"ctgA");
        load(&decoder, b"ctgB");
        decoder.retain(&HashSet::from([&b"ctgA"[..]])).unwrap();
        assert_eq!(held(&decoder), [b"ctgA".to_vec(), b"ctgB".to_vec()]);
        decoder.retain(&HashSet::from([&b"ctgB"[..]])).unwrap();
        assert_eq!(held(&decoder), [b"ctgB".to_vec()]);
        assert_eq!(decoder.repository.len(), 1);
        assert_eq!(decoder.held_sequence(b"ctgB").unwrap().len(), 200);
        decoder.retain(&HashSet::new()).unwrap();
        assert!(held(&decoder).is_empty());
        assert!(decoder.repository.is_empty());
        assert!(decoder.held_sequence(b"ctgA").is_none());
        assert!(decoder.repository.is_empty());
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
    fn test_md_nm_stop_at_the_end_of_the_reference_like_htslib() {
        assert_eq!(
            regenerated("r\t0\tc\t27\t60\t6M\t*\t0\t0\tACGTAA\t*"),
            "MD:Z:4 NM:i:0"
        );
        assert_eq!(
            regenerated("r\t0\tc\t27\t60\t6M\t*\t0\t0\tAGGTAA\t*"),
            "MD:Z:1C2 NM:i:1"
        );
        assert_eq!(
            regenerated("r\t0\tc\t29\t60\t2M2I3M\t*\t0\t0\tGTAAACG\t*"),
            "MD:Z:2 NM:i:2"
        );
        assert_eq!(
            regenerated("r\t0\tc\t28\t60\t2M3D1M\t*\t0\t0\tCGA\t*"),
            "MD:Z:2^T0 NM:i:1"
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
    fn test_a_reference_shorter_than_the_header_length_fails_cleanly() {
        let path = test_data("contigs.cram");
        let mut reader = cram::io::Reader::new(File::open(&path).unwrap());
        let header = reader.read_header().unwrap();
        let reference = Reference::load(&test_data("contigs.short.fa")).unwrap();
        let mut decoder = CramDecoder::new(Some(reference), &path, &header);
        let mut container = cram::io::reader::Container::default();
        reader.read_container(&mut container).unwrap();
        let error = decoder.decode(&container, &header).unwrap_err();
        assert!(format!("{error:#}").contains(
            "reference sequence `ctgA` is 300 bases in the --reference FASTA but 400 in the @SQ \
             header"
        ));
    }

    #[test]
    fn test_reference_reads_a_whole_sequence() {
        let reference = Reference::load(&test_data("contigs.fa")).unwrap();
        assert_eq!(reference.sequence(b"ctgB").unwrap().len(), 200);
        assert!(reference.sequence(b"missing").is_err());
    }
}
