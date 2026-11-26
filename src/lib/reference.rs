//! CRAM decoding against a `--reference` FASTA, holding only the reference sequences slices need.

use std::collections::{HashMap, HashSet};
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, RwLock, RwLockReadGuard};

use anyhow::{anyhow, bail, Context, Result};
use noodles::core::Region;
use noodles::fasta::fai;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::{Record as _, RecordBuf};
use noodles::{bgzf, cram, fasta, sam};

use crate::cram_raw::{self, EmbeddedReference, Recorder};
use crate::reference_check::HeaderCheck;

/// An indexed FASTA (with a `.fai`) that CRAM inputs decode against; clones share one reader.
#[derive(Clone)]
pub struct Reference {
    reader: Arc<Mutex<fasta::io::IndexedReader<fasta::io::BufReader<File>>>>,
    check: Option<Arc<HeaderCheck>>,
}

impl Reference {
    /// Open an indexed FASTA (`PATH.fai` alongside); `check` enables [`Self::check_header`].
    pub fn open(path: &Path, check: bool, threads: usize) -> Result<Self> {
        let context = || format!("failed to open the --reference FASTA {}", path.display());
        File::open(path).with_context(context)?;
        let mut index_path = OsString::from(path);
        index_path.push(".fai");
        let index_path = PathBuf::from(index_path);
        let index = fai::fs::read(&index_path).with_context(|| {
            format!(
                "failed to read the .fai index {} of the --reference FASTA {}; create it with samtools faidx",
                index_path.display(),
                path.display()
            )
        })?;
        reject_duplicate_names(&index, &index_path)?;
        if let Some("gz" | "bgz") = path.extension().and_then(|ext| ext.to_str()) {
            let mut gzi_path = OsString::from(path);
            gzi_path.push(".gzi");
            let gzi_path = PathBuf::from(gzi_path);
            bgzf::gzi::fs::read(&gzi_path).with_context(|| {
                format!(
                    "failed to read the .gzi index {} of the bgzipped --reference FASTA {}; create it with bgzip -r",
                    gzi_path.display(),
                    path.display()
                )
            })?;
        }
        let reader = fasta::io::indexed_reader::Builder::default()
            .set_index(index.clone())
            .build_from_path(path)
            .with_context(context)?;
        Ok(Self {
            reader: Arc::new(Mutex::new(reader)),
            check: check.then(|| Arc::new(HeaderCheck::new(path, index, threads))),
        })
    }

    /// Fail unless every `@SQ` contig is in this FASTA with the same `LN` and any `M5`.
    pub fn check_header(&self, header: &sam::Header, input: &Path) -> Result<()> {
        match &self.check {
            Some(check) => check.check(header, input),
            None => Ok(()),
        }
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

/// Fail when the `.fai` lists a name twice, since a lookup by name could pick either record.
fn reject_duplicate_names(index: &fai::Index, path: &Path) -> Result<()> {
    let mut seen = HashSet::new();
    for record in index.as_ref() {
        let name: &[u8] = record.name().as_ref();
        if !seen.insert(name) {
            bail!(
                "contig `{}` is listed more than once in the .fai index {}; give each sequence of the --reference FASTA a unique name",
                String::from_utf8_lossy(name),
                path.display()
            );
        }
    }
    Ok(())
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

/// Serves zeroed stand-ins that noodles indexes by position in a slice that needs no reference.
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

/// Decodes one CRAM input's containers into owned records. Clones share one set of held sequences,
/// so decode workers never copy a sequence another worker already holds.
#[derive(Clone)]
pub(crate) struct CramDecoder {
    shared: Arc<Shared>,
}

/// The state [`CramDecoder`] clones share. Slices decode under a read lock on the repository, and
/// dropping held sequences takes the write lock, so it waits for the slices in flight.
struct Shared {
    reference: Option<Reference>,
    path: PathBuf,
    lengths: Arc<HashMap<Vec<u8>, usize>>,
    repository: RwLock<fasta::Repository>,
    held: Held,
    pinned: Mutex<HashMap<Vec<u8>, usize>>,
    newest: Mutex<HashMap<Vec<u8>, usize>>,
    tickets: AtomicUsize,
    blank: Mutex<Blank>,
}

/// The zeroed stand-ins slices written without the reference decode against, and their lengths.
struct Blank {
    repository: Option<fasta::Repository>,
    lengths: HashMap<Vec<u8>, usize>,
}

/// Marks a sequence as wanted by a slice in flight until dropped, so other workers keep it.
struct Pin<'a> {
    pinned: &'a Mutex<HashMap<Vec<u8>, usize>>,
    name: Vec<u8>,
}

impl Drop for Pin<'_> {
    fn drop(&mut self) {
        if let Ok(mut pinned) = self.pinned.lock() {
            if let Some(count) = pinned.get_mut(&self.name) {
                *count -= 1;
                if *count == 0 {
                    pinned.remove(&self.name);
                }
            }
        }
    }
}

impl CramDecoder {
    /// A decoder for the CRAM at `path` with `header`, reading `reference` sequences when given.
    pub(crate) fn new(reference: Option<Reference>, path: &Path, header: &sam::Header) -> Self {
        let lengths: HashMap<Vec<u8>, usize> = header
            .reference_sequences()
            .iter()
            .map(|(name, sequence)| (name.to_vec(), sequence.length().get()))
            .collect();
        let lengths = Arc::new(lengths);
        let held = Held::default();
        let repository = loading_repository(&reference, path, &lengths, &held, HashMap::new());
        Self {
            shared: Arc::new(Shared {
                reference,
                path: path.to_path_buf(),
                lengths: Arc::clone(&lengths),
                repository: RwLock::new(repository),
                held,
                pinned: Mutex::default(),
                newest: Mutex::default(),
                tickets: AtomicUsize::new(0),
                blank: Mutex::new(Blank {
                    repository: None,
                    lengths: (*lengths).clone(),
                }),
            }),
        }
    }

    fn held(&self) -> Result<std::sync::MutexGuard<'_, HashMap<Vec<u8>, usize>>> {
        self.shared
            .held
            .lock()
            .map_err(|_| anyhow!("the CRAM reference bookkeeping was poisoned"))
    }

    fn pinned(&self) -> Result<std::sync::MutexGuard<'_, HashMap<Vec<u8>, usize>>> {
        self.shared
            .pinned
            .lock()
            .map_err(|_| anyhow!("the CRAM reference bookkeeping was poisoned"))
    }

    fn blank(&self) -> Result<std::sync::MutexGuard<'_, Blank>> {
        self.shared.blank.lock().map_err(|_| poisoned())
    }

    /// Record that container `order` uses `names`, so no earlier container drops them.
    fn mark<'a>(&self, names: impl IntoIterator<Item = &'a [u8]>, order: usize) -> Result<()> {
        let mut newest = self.shared.newest.lock().map_err(|_| poisoned())?;
        for name in names {
            let last = newest.entry(name.to_vec()).or_default();
            *last = (*last).max(order);
        }
        Ok(())
    }

    fn pin(&self, name: &[u8]) -> Result<Pin<'_>> {
        *self.pinned()?.entry(name.to_vec()).or_default() += 1;
        Ok(Pin {
            pinned: &self.shared.pinned,
            name: name.to_vec(),
        })
    }

    fn read(&self) -> Result<RwLockReadGuard<'_, fasta::Repository>> {
        self.shared.repository.read().map_err(|_| poisoned())
    }

    /// A held sequence, never loading one.
    fn held_sequence(
        &self,
        repository: &fasta::Repository,
        name: &[u8],
    ) -> Option<Arc<fasta::record::Sequence>> {
        if !self.held().ok()?.contains_key(name) {
            return None;
        }
        repository.get(name).and_then(Result::ok)
    }

    /// Decode every slice of `container` into owned records, in file order.
    ///
    /// Each slice decodes once against every sequence it references. A probe with no sequences
    /// first finds the first one it needs, so a single-reference container, or sorted input, drops
    /// what the slice cannot use before loading more; afterwards the unused sequences are dropped.
    /// A sequence another clone's slice still needs, or a later container already used, is kept. A
    /// container written without the reference never loads one, and a slice that embeds its
    /// reference regenerates MD and NM from it, as htslib does.
    pub(crate) fn decode(
        &self,
        container: &cram::io::reader::Container,
        header: &sam::Header,
        raw: &[u8],
    ) -> Result<Vec<RecordBuf>> {
        self.decode_in_order(self.ticket(), container, header, raw)
    }

    /// The place in file order of the next container read, for [`Self::decode_in_order`].
    pub(crate) fn ticket(&self) -> usize {
        self.shared.tickets.fetch_add(1, Ordering::Relaxed)
    }

    /// Decode `container`, the one at place `order` in file order, as [`Self::decode`] does.
    pub(crate) fn decode_in_order(
        &self,
        order: usize,
        container: &cram::io::reader::Container,
        header: &sam::Header,
        raw: &[u8],
    ) -> Result<Vec<RecordBuf>> {
        quietly(|| self.decode_container(order, container, header, raw)).unwrap_or_else(|| {
            Err(anyhow!(
                "failed to decode a CRAM container of {}: it is malformed or uses a CRAM feature \
                 unmux cannot decode",
                self.shared.path.display()
            ))
        })
    }

    fn decode_container(
        &self,
        order: usize,
        container: &cram::io::reader::Container,
        header: &sam::Header,
        raw: &[u8],
    ) -> Result<Vec<RecordBuf>> {
        let compression_header = container
            .compression_header()
            .context("failed to read CRAM compression header")?;
        let regenerate = uses_external_reference(&compression_header)?;
        let single_reference = !container.header().reference_sequence_context().is_many();
        let sorted = is_coordinate_sorted(header);
        let mut out = Vec::new();
        let landmarks = container.header().landmarks();
        for (index, slice) in container.slices().enumerate() {
            let slice = slice.context("failed to read CRAM slice")?;
            let (core, external) = slice
                .decode_blocks()
                .context("failed to decode CRAM slice blocks")?;
            let embedded = embedded_bases(raw, landmarks, index, &external)?;
            let decode = |repository: fasta::Repository| {
                slice.records(repository, header, &compression_header, &core, &external)
            };
            let mut pin = None;
            let mut probed = None;
            if regenerate {
                let result = decode(fasta::Repository::new(ProbeAdapter));
                if let Some(name) = result.as_ref().err().and_then(requested) {
                    self.mark([name], order)?;
                    pin = Some(self.pin(name)?);
                    if single_reference {
                        self.retain(&HashSet::from([name]), order)?;
                    } else if sorted {
                        self.retain_from(header, name, order)?;
                    }
                }
                probed = Some(result);
            }
            let repository = self.read()?;
            let records = match probed {
                Some(Err(error)) if requested(&error).is_some() => decode(repository.clone()),
                Some(probed) => probed,
                None => {
                    let records = decode(self.blank_repository(&[])?)
                        .context("failed to read CRAM slice records")?;
                    let ends: Vec<(&[u8], usize)> = records
                        .iter()
                        .map(|record| aligned_end(header, record))
                        .filter_map(Result::transpose)
                        .collect::<Result<_>>()?;
                    if sorted {
                        self.blank()?.repository = None;
                    }
                    if let Some((embedded, bases)) = embedded {
                        let end = ends.iter().map(|&(_, end)| end).max().unwrap_or(0);
                        if end >= embedded.start + bases.len() {
                            drop(records);
                            out.extend(self.decode_past_embedded(
                                (raw, landmarks),
                                (header, &compression_header),
                                index,
                                (embedded, bases),
                                end,
                            )?);
                            continue;
                        }
                    }
                    let covered = {
                        let blank = self.blank()?;
                        ends.iter()
                            .all(|(name, end)| blank.lengths.get(*name) >= Some(end))
                    };
                    if covered {
                        Ok(records)
                    } else {
                        drop(records);
                        decode(self.blank_repository(&ends)?)
                    }
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
                decode(self.padded_repository(&repository, &padded)?)
                    .context("failed to read CRAM slice records")?
            };
            let used = records
                .iter()
                .filter_map(|record| sequence_name(header, record).transpose())
                .collect::<Result<HashSet<_>>>()?;
            self.mark(used.iter().copied(), order)?;
            for (i, record) in records.iter().enumerate() {
                let name = sequence_name(header, record)?;
                let record = match name.filter(|name| padded.contains_key(name)) {
                    Some(_) => &repadded[i],
                    None => record,
                };
                let held = name
                    .filter(|_| regenerate && embedded.is_none())
                    .and_then(|name| self.held_sequence(&repository, name));
                let reference = match (embedded, &held) {
                    (Some((embedded, bases)), _) => Some((bases, embedded.start)),
                    (None, Some(sequence)) => Some(((**sequence).as_ref(), 1)),
                    (None, None) => None,
                };
                let converted = if regenerate {
                    owned(header, record, reference)
                } else {
                    quietly(|| owned(header, record, reference))
                        .unwrap_or_else(|| Err(self.unreferenced(record)))
                }?;
                // A zeroed stand-in base means the record needed reference bases its slice lacks.
                if !regenerate && converted.sequence().as_ref().contains(&0) {
                    return Err(self.unreferenced(record));
                }
                out.push(converted);
            }
            drop((records, repadded, repository, pin));
            self.retain(&used, order)?;
        }
        Ok(out)
    }

    /// Decode slice `index` of container `raw` with its embedded reference padded with N to `end`.
    fn decode_past_embedded(
        &self,
        (raw, landmarks): (&[u8], &[usize]),
        (header, compression_header): (&sam::Header, &cram::container::CompressionHeader),
        index: usize,
        (embedded, bases): (EmbeddedReference, &[u8]),
        end: usize,
    ) -> Result<Vec<RecordBuf>> {
        let mut padded = bases.to_vec();
        padded.resize(end + 1 - embedded.start, b'N');
        let rebuilt = cram_raw::with_embedded_bases(raw, landmarks, index, embedded, &padded)
            .context("failed to lengthen an embedded CRAM reference")?;
        let mut container = cram::io::reader::Container::default();
        cram::io::Reader::new(&rebuilt[..])
            .read_container(&mut container)
            .context("failed to lengthen an embedded CRAM reference")?;
        let slice = container
            .slices()
            .nth(index)
            .ok_or_else(|| anyhow!("failed to lengthen an embedded CRAM reference"))?
            .context("failed to read CRAM slice")?;
        let (core, external) = slice
            .decode_blocks()
            .context("failed to decode CRAM slice blocks")?;
        let repository = fasta::Repository::default();
        let records = slice
            .records(repository, header, compression_header, &core, &external)
            .context("failed to read CRAM slice records")?;
        records
            .iter()
            .map(|record| {
                quietly(|| owned(header, record, Some((bases, embedded.start))))
                    .unwrap_or_else(|| Err(self.unreferenced(record)))
            })
            .collect()
    }

    /// The error for a record that needs reference bases its slice does not carry.
    fn unreferenced(&self, record: &cram::Record<'_>) -> anyhow::Error {
        anyhow!(
            "CRAM record `{}` in {} needs reference bases its slice does not carry",
            record_name(record),
            self.shared.path.display()
        )
    }

    /// Zeroed stand-ins as long as each header length, lengthened to reach every one of `ends`.
    fn blank_repository(&self, ends: &[(&[u8], usize)]) -> Result<fasta::Repository> {
        let mut blank = self.blank()?;
        for &(name, end) in ends {
            let length = blank.lengths.entry(name.to_vec()).or_default();
            if end > *length {
                *length = end;
                blank.repository = None;
            }
        }
        let lengths = blank.lengths.clone();
        Ok(blank
            .repository
            .get_or_insert_with(|| fasta::Repository::new(BlankAdapter(lengths)))
            .clone())
    }

    /// Held sequences extended with N to the given lengths, and an empty stand-in for any other.
    fn padded_repository(
        &self,
        repository: &fasta::Repository,
        padded: &HashMap<&[u8], usize>,
    ) -> Result<fasta::Repository> {
        let mut sequences = HashMap::new();
        for (&name, &length) in padded {
            let held = self
                .held_sequence(repository, name)
                .ok_or_else(|| anyhow!("a padded CRAM reference sequence was not held"))?;
            let mut bases = Vec::with_capacity(length);
            bases.extend_from_slice((*held).as_ref());
            bases.resize(length.max(bases.len()), b'N');
            sequences.insert(name.to_vec(), fasta::record::Sequence::from(bases));
        }
        Ok(fasta::Repository::new(PaddedAdapter(sequences)))
    }

    /// Drop the held sequences ahead of `first` in header order, which sorted input never revisits.
    fn retain_from(&self, header: &sam::Header, first: &[u8], order: usize) -> Result<()> {
        let sequences = header.reference_sequences();
        let first = sequences.get_index_of(first);
        let needed: Vec<Vec<u8>> = self
            .held()?
            .keys()
            .filter(|name| sequences.get_index_of(name.as_slice()) >= first)
            .cloned()
            .collect();
        self.retain(&needed.iter().map(Vec::as_slice).collect(), order)
    }

    /// The held sequences container `order` drops to keep only `needed`, what slices in flight
    /// pinned and what later containers used, or none when they free less than keeping the rest
    /// copies.
    fn droppable(&self, needed: &HashSet<&[u8]>, order: usize) -> Result<Vec<Vec<u8>>> {
        let held = self.held()?.clone();
        let pinned = self.pinned()?;
        let newest = self.shared.newest.lock().map_err(|_| poisoned())?;
        let keep = |name: &Vec<u8>| {
            needed.contains(name.as_slice())
                || pinned.contains_key(name)
                || newest.get(name).is_some_and(|&last| last > order)
        };
        let (kept, dropped): (Vec<_>, Vec<_>) = held.iter().partition(|(name, _)| keep(name));
        let bases = |part: &[(&Vec<u8>, &usize)]| part.iter().map(|(_, len)| **len).sum::<usize>();
        if dropped.is_empty() || bases(&dropped) < bases(&kept) {
            return Ok(Vec::new());
        }
        Ok(dropped.into_iter().map(|(name, _)| name.clone()).collect())
    }

    /// Drop the held sequences container `order` does not need, as [`Self::droppable`] picks them.
    fn retain(&self, needed: &HashSet<&[u8]>, order: usize) -> Result<()> {
        if self.droppable(needed, order)?.is_empty() {
            return Ok(());
        }
        let mut repository = self.shared.repository.write().map_err(|_| poisoned())?;
        let dropped = self.droppable(needed, order)?;
        if dropped.is_empty() {
            return Ok(());
        }
        let kept: Vec<Vec<u8>> = self
            .held()?
            .keys()
            .filter(|name| !dropped.contains(name))
            .cloned()
            .collect();
        let carried: HashMap<Vec<u8>, Arc<fasta::record::Sequence>> = kept
            .into_iter()
            .filter_map(|name| Some((name.clone(), repository.get(&name)?.ok()?)))
            .collect();
        let names: Vec<Vec<u8>> = carried.keys().cloned().collect();
        self.held()?.clear();
        *repository = fasta::Repository::default();
        *repository = loading_repository(
            &self.shared.reference,
            &self.shared.path,
            &self.shared.lengths,
            &self.shared.held,
            carried,
        );
        for name in names {
            repository
                .get(&name)
                .transpose()
                .context("failed to keep a CRAM reference sequence")?;
        }
        Ok(())
    }
}

fn loading_repository(
    reference: &Option<Reference>,
    path: &Path,
    lengths: &Arc<HashMap<Vec<u8>, usize>>,
    held: &Held,
    carried: HashMap<Vec<u8>, Arc<fasta::record::Sequence>>,
) -> fasta::Repository {
    fasta::Repository::new(LoadingAdapter {
        reference: reference.clone(),
        path: path.to_path_buf(),
        lengths: Arc::clone(lengths),
        carried,
        held: Arc::clone(held),
    })
}

fn poisoned() -> anyhow::Error {
    anyhow!("a CRAM decode thread panicked")
}

/// Read the next container of `reader` into `container`, keeping its raw bytes in the reader's
/// [`Recorder`]; `false` at the end of the input.
pub(crate) fn read_container<R: io::Read>(
    reader: &mut cram::io::Reader<Recorder<R>>,
    container: &mut cram::io::reader::Container,
) -> io::Result<bool> {
    reader.get_mut().clear();
    Ok(reader.read_container(container)? != 0)
}

/// The reference the slice at `index` of container `raw` embeds, with its bases from `external`.
fn embedded_bases<'a>(
    raw: &[u8],
    landmarks: &[usize],
    index: usize,
    external: &'a [(i32, std::borrow::Cow<'a, [u8]>)],
) -> Result<Option<(EmbeddedReference, &'a [u8])>> {
    let slice = cram_raw::slice(raw, landmarks, index).context("failed to read a CRAM slice")?;
    let Some(embedded) =
        cram_raw::embedded_reference(slice).context("failed to read a CRAM slice header")?
    else {
        return Ok(None);
    };
    Ok(external
        .iter()
        .find(|(id, _)| *id == embedded.content_id)
        .map(|(_, bases)| (embedded, bases.as_ref())))
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

/// A record's name, for messages.
fn record_name(record: &cram::Record<'_>) -> String {
    record.name().map(ToString::to_string).unwrap_or_default()
}

thread_local! {
    static QUIET: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Run `f`, returning `None` without printing the panic if noodles panics on input it cannot decode.
pub(crate) fn quietly<T>(f: impl FnOnce() -> T) -> Option<T> {
    static HOOK: std::sync::Once = std::sync::Once::new();
    HOOK.call_once(|| {
        let previous = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            if !QUIET.with(std::cell::Cell::get) {
                previous(info);
            }
        }));
    });
    let quiet = QUIET.with(|quiet| quiet.replace(true));
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(f));
    QUIET.with(|cell| cell.set(quiet));
    result.ok()
}

/// The reference bases a record decoded against and the position of the first.
type Bases<'a> = Option<(&'a [u8], usize)>;

/// Convert a decoded CRAM record into an owned record, as htslib would decode it.
fn owned(
    header: &sam::Header,
    record: &cram::Record<'_>,
    reference: Bases<'_>,
) -> Result<RecordBuf> {
    let mut record = RecordBuf::try_from_alignment_record(header, record)
        .context("failed to convert CRAM record")?;
    // CRAM stores no MAPQ for an unmapped record; noodles reports 255, htslib 0.
    if record.flags().is_unmapped() && record.mapping_quality().is_none() {
        *record.mapping_quality_mut() = MappingQuality::new(0);
    }
    record.sequence_mut().as_mut().make_ascii_uppercase();
    regenerate_md_nm(header, &mut record, reference);
    Ok(record)
}

/// Add the MD and NM tags htslib regenerates when it decodes a mapped record against `reference`.
fn regenerate_md_nm(header: &sam::Header, record: &mut RecordBuf, reference: Bases<'_>) {
    let skip = take_cram_flags(record);
    let want_md = skip & 1 == 0 && record.data().get(&Tag::MISMATCHED_POSITIONS).is_none();
    let want_nm = skip & 2 == 0 && record.data().get(&Tag::EDIT_DISTANCE).is_none();
    if !(want_md || want_nm) || record.flags().is_unmapped() || record.sequence().is_empty() {
        return;
    }
    let (Some((bases, start)), Some(length)) = (
        reference,
        record
            .reference_sequence_id()
            .and_then(|id| header.reference_sequences().get_index(id))
            .map(|(_, map)| map.length().get()),
    ) else {
        return;
    };
    let within = length.saturating_sub(start - 1).min(bases.len());
    if let Some((md, nm)) = md_nm(record, &bases[..within], start - 1) {
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

/// The MD string and edit distance htslib derives from a record's CIGAR, SEQ and the reference
/// bases after `offset`; like htslib, the MD string stops at the first base past their end.
fn md_nm(record: &RecordBuf, reference: &[u8], offset: usize) -> Option<(Vec<u8>, u32)> {
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
                    let base = ref_pos.checked_sub(offset).and_then(|i| reference.get(i));
                    match base.map(u8::to_ascii_uppercase) {
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
                    let at = ref_pos.checked_sub(offset);
                    if let Some(bases) =
                        at.and_then(|at| reference.get(at..reference.len().min(at + len)))
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
        decode_with(
            name,
            Some(Reference::open(&test_data("contigs.fa"), false, 1).unwrap()),
        )
    }

    /// Decode a test CRAM against `reference`, as [`decode_all`] does.
    fn decode_with(name: &str, reference: Option<Reference>) -> (Vec<RecordBuf>, usize) {
        decode_file(&test_data(name), reference)
    }

    /// Decode the CRAM at `path` against `reference`, as [`decode_all`] does.
    fn decode_file(path: &Path, reference: Option<Reference>) -> (Vec<RecordBuf>, usize) {
        let mut reader = cram::io::Reader::new(Recorder::new(File::open(path).unwrap()));
        let header = reader.read_header().unwrap();
        let decoder = CramDecoder::new(reference, path, &header);
        let mut container = cram::io::reader::Container::default();
        let mut records = Vec::new();
        let mut most = 0;
        while read_container(&mut reader, &mut container).unwrap() {
            records.extend(
                decoder
                    .decode(&container, &header, reader.get_ref().bytes())
                    .unwrap(),
            );
            most = most.max(decoder.read().unwrap().len());
        }
        (records, most)
    }

    #[test]
    fn test_the_name_star_decodes_as_a_missing_name() {
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("star.cram");
        let header = sam::Header::default();
        let mut writer = cram::io::Writer::new(File::create(&path).unwrap());
        writer.write_header(&header).unwrap();
        for name in ["n1", "*", "n3"] {
            let record = RecordBuf::builder()
                .set_name(name)
                .set_flags(Flags::UNMAPPED)
                .set_sequence(Sequence::from(b"ACGT".to_vec()))
                .set_quality_scores(QualityScores::from(vec![30; 4]))
                .build();
            writer.write_alignment_record(&header, &record).unwrap();
        }
        writer.try_finish(&header).unwrap();

        let (records, _) = decode_file(&path, None);
        let names: Vec<Option<String>> = records
            .iter()
            .map(|record| record.name().map(ToString::to_string))
            .collect();
        assert_eq!(
            names,
            [Some("n1".to_string()), None, Some("n3".to_string())]
        );
    }

    #[test]
    fn test_the_reference_flag_is_read_from_the_compression_header() {
        for (name, external) in [
            ("contigs.cram", true),
            ("contigs.multi-ref.cram", true),
            ("contigs.no-ref.cram", false),
            ("contigs.embedded-ref.cram", false),
        ] {
            let mut reader =
                cram::io::Reader::new(Recorder::new(File::open(test_data(name)).unwrap()));
            reader.read_header().unwrap();
            let mut container = cram::io::reader::Container::default();
            read_container(&mut reader, &mut container).unwrap();
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
            Some(Reference::open(&test_data("contigs.fa"), false, 1).unwrap()),
        ] {
            let mut reader = cram::io::Reader::new(Recorder::new(File::open(&path).unwrap()));
            let header = reader.read_header().unwrap();
            let decoder = CramDecoder::new(reference, &path, &header);
            let mut container = cram::io::reader::Container::default();
            let mut records = Vec::new();
            while read_container(&mut reader, &mut container).unwrap() {
                records.extend(
                    decoder
                        .decode(&container, &header, reader.get_ref().bytes())
                        .unwrap(),
                );
            }
            assert_eq!(records.len(), 29);
            assert!(decoder.read().unwrap().is_empty());
            assert!(records
                .iter()
                .all(|record| !record.sequence().as_ref().contains(&0)));
        }
    }

    fn tagged(records: &[RecordBuf], tag: Tag) -> usize {
        records
            .iter()
            .filter(|record| record.data().get(&tag).is_some())
            .count()
    }

    #[test]
    fn test_embedded_reference_slices_get_md_and_nm_like_htslib() {
        let (embedded, _) = decode_with("contigs.embedded-ref.cram", None);
        let (external, _) = decode_all("contigs.cram");
        assert_eq!(tagged(&embedded, Tag::MISMATCHED_POSITIONS), 26);
        assert_eq!(embedded, external);
    }

    #[test]
    fn test_embedded_consensus_reference_slices_keep_md_and_nm_absent() {
        let (records, _) = decode_with("contigs.embedded-ref-2.cram", None);
        assert_eq!(records.len(), 29);
        assert_eq!(tagged(&records, Tag::MISMATCHED_POSITIONS), 1);
        assert_eq!(tagged(&records, Tag::new(b'c', b'F')), 0);
    }

    #[test]
    fn test_records_past_the_end_of_an_embedded_reference_decode_like_htslib() {
        let (embedded, _) = decode_with("past-end.embedded-ref.cram", None);
        let (external, _) = decode_all("past-end.cram");
        assert_eq!(embedded.len(), 5);
        assert_eq!(embedded, external);
    }

    #[test]
    fn test_embedded_reference_is_read_from_the_slice_header() {
        let path = test_data("contigs.embedded-ref.slices.cram");
        let mut reader = cram::io::Reader::new(Recorder::new(File::open(path).unwrap()));
        reader.read_header().unwrap();
        let mut container = cram::io::reader::Container::default();
        read_container(&mut reader, &mut container).unwrap();
        let raw = reader.get_ref().bytes();
        let landmarks = container.header().landmarks();
        let starts: Vec<usize> = (0..landmarks.len())
            .map(|i| {
                let slice = cram_raw::slice(raw, landmarks, i).unwrap();
                cram_raw::embedded_reference(slice).unwrap().unwrap().start
            })
            .collect();
        assert_eq!(starts, [11, 51, 91, 121]);
    }

    #[test]
    fn test_a_container_rebuilt_with_longer_embedded_bases_keeps_every_slice() {
        let path = test_data("contigs.embedded-ref.slices.cram");
        let mut reader = cram::io::Reader::new(Recorder::new(File::open(path).unwrap()));
        let header = reader.read_header().unwrap();
        let mut container = cram::io::reader::Container::default();
        read_container(&mut reader, &mut container).unwrap();
        let raw = reader.get_ref().bytes().to_vec();
        let landmarks = container.header().landmarks().to_vec();
        let decode = |container: &cram::io::reader::Container| {
            let compression_header = container.compression_header().unwrap();
            let mut names = Vec::new();
            let mut embedded = Vec::new();
            for slice in container.slices() {
                let slice = slice.unwrap();
                let (core, external) = slice.decode_blocks().unwrap();
                let records = slice
                    .records(
                        fasta::Repository::default(),
                        &header,
                        &compression_header,
                        &core,
                        &external,
                    )
                    .unwrap();
                names.extend(records.iter().map(record_name));
                embedded.push(
                    external
                        .iter()
                        .find(|(id, _)| *id == 10)
                        .unwrap()
                        .1
                        .to_vec(),
                );
            }
            (names, embedded)
        };
        let (names, mut embedded) = decode(&container);
        let slice = cram_raw::slice(&raw, &landmarks, 1).unwrap();
        let reference = cram_raw::embedded_reference(slice).unwrap().unwrap();
        embedded[1].extend_from_slice(b"NNNNNNN");
        let rebuilt =
            cram_raw::with_embedded_bases(&raw, &landmarks, 1, reference, &embedded[1]).unwrap();
        let mut container = cram::io::reader::Container::default();
        cram::io::Reader::new(&rebuilt[..])
            .read_container(&mut container)
            .unwrap();
        assert_eq!(decode(&container), (names, embedded));
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
        let reference = Reference::open(&test_data("contigs.fa"), false, 1).unwrap();
        let path = test_data("contigs.cram");
        let header = cram::io::Reader::new(Recorder::new(File::open(&path).unwrap()))
            .read_header()
            .unwrap();
        let decoder = CramDecoder::new(Some(reference), &path, &header);
        let load = |decoder: &CramDecoder, name: &[u8]| {
            decoder.read().unwrap().get(name).unwrap().unwrap();
        };
        let held = |decoder: &CramDecoder| {
            let mut names: Vec<Vec<u8>> = decoder.held().unwrap().keys().cloned().collect();
            names.sort();
            names
        };
        load(&decoder, b"ctgA");
        load(&decoder, b"ctgB");
        decoder.retain(&HashSet::from([&b"ctgA"[..]]), 0).unwrap();
        assert_eq!(held(&decoder), [b"ctgA".to_vec(), b"ctgB".to_vec()]);
        decoder.retain(&HashSet::from([&b"ctgB"[..]]), 0).unwrap();
        assert_eq!(held(&decoder), [b"ctgB".to_vec()]);
        let repository = decoder.read().unwrap().clone();
        assert_eq!(repository.len(), 1);
        assert_eq!(
            decoder.held_sequence(&repository, b"ctgB").unwrap().len(),
            200
        );
        decoder.retain(&HashSet::new(), 0).unwrap();
        let repository = decoder.read().unwrap().clone();
        assert!(held(&decoder).is_empty());
        assert!(repository.is_empty());
        assert!(decoder.held_sequence(&repository, b"ctgA").is_none());
        assert!(repository.is_empty());
    }

    /// A decoder for the test CRAM `name` against the test FASTA.
    fn test_decoder(name: &str) -> CramDecoder {
        let path = test_data(name);
        let header = cram::io::Reader::new(Recorder::new(File::open(&path).unwrap()))
            .read_header()
            .unwrap();
        let reference = Reference::open(&test_data("contigs.fa"), false, 1).unwrap();
        CramDecoder::new(Some(reference), &path, &header)
    }

    #[test]
    fn test_retain_keeps_sequences_pinned_by_slices_in_flight() {
        let decoder = test_decoder("contigs.cram");
        for name in [&b"ctgA"[..], b"ctgB"] {
            decoder.read().unwrap().get(name).unwrap().unwrap();
        }
        let pin = decoder.pin(b"ctgA").unwrap();
        decoder.retain(&HashSet::new(), 0).unwrap();
        assert_eq!(decoder.held().unwrap().len(), 2);
        drop(pin);
        decoder.retain(&HashSet::new(), 0).unwrap();
        assert!(decoder.held().unwrap().is_empty());
    }

    #[test]
    fn test_retain_keeps_sequences_a_later_container_used() {
        let decoder = test_decoder("contigs.cram");
        for name in [&b"ctgA"[..], b"ctgB"] {
            decoder.read().unwrap().get(name).unwrap().unwrap();
        }
        decoder.mark([&b"ctgB"[..]], 5).unwrap();
        decoder.retain(&HashSet::new(), 3).unwrap();
        let held: Vec<Vec<u8>> = decoder.held().unwrap().keys().cloned().collect();
        assert_eq!(held, [b"ctgB".to_vec()]);
        decoder.retain(&HashSet::new(), 5).unwrap();
        assert!(decoder.held().unwrap().is_empty());
    }

    #[test]
    fn test_decoder_clones_decode_concurrently_like_one_decoder() {
        for name in [
            "contigs.cram",
            "contigs.multi-ref.cram",
            "contigs.mixed-ref.cram",
            "past-end.cram",
            "past-end.embedded-ref.cram",
            "contigs.embedded-ref.slices.cram",
            "contigs.no-ref.multi-ref.cram",
        ] {
            let (single, _) = decode_all(name);
            let mut reader =
                cram::io::Reader::new(Recorder::new(File::open(test_data(name)).unwrap()));
            let header = reader.read_header().unwrap();
            let mut containers = Vec::new();
            let mut container = cram::io::reader::Container::default();
            while read_container(&mut reader, &mut container).unwrap() {
                containers.push((std::mem::take(&mut container), reader.get_mut().take()));
            }
            let decoder = test_decoder(name);
            let decoded: Vec<Vec<RecordBuf>> = std::thread::scope(|scope| {
                let handles: Vec<_> = containers
                    .iter()
                    .map(|(container, raw)| {
                        let decoder = decoder.clone();
                        let header = &header;
                        scope.spawn(move || decoder.decode(container, header, raw).unwrap())
                    })
                    .collect();
                handles.into_iter().map(|h| h.join().unwrap()).collect()
            });
            assert_eq!(decoded.concat(), single, "{name}");
            assert!(decoder.held().unwrap().len() <= 3, "{name}");
            assert!(decoder.pinned().unwrap().is_empty(), "{name}");
        }
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
        regenerate_md_nm(&header, &mut record, Some((CONTIG, 1)));
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
    fn test_md_nm_read_reference_bases_that_start_past_position_one() {
        let (header, mut record) = parse("r\t0\tc\t9\t60\t6M\t*\t0\t0\tGTACTA\t*");
        regenerate_md_nm(&header, &mut record, Some((&CONTIG[6..], 7)));
        assert_eq!(
            record.data().get(&Tag::MISMATCHED_POSITIONS),
            Some(&Value::String(b"4G0T0".into()))
        );
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
        let mut reader = cram::io::Reader::new(Recorder::new(File::open(&path).unwrap()));
        let header = reader.read_header().unwrap();
        let reference = Reference::open(&test_data("contigs.short.fa"), false, 1).unwrap();
        let decoder = CramDecoder::new(Some(reference), &path, &header);
        let mut container = cram::io::reader::Container::default();
        read_container(&mut reader, &mut container).unwrap();
        let error = decoder
            .decode(&container, &header, reader.get_ref().bytes())
            .unwrap_err();
        assert!(format!("{error:#}").contains(
            "reference sequence `ctgA` is 300 bases in the --reference FASTA but 400 in the @SQ \
             header"
        ));
    }

    #[test]
    fn test_a_panic_inside_noodles_becomes_an_error() {
        let path = test_data("contigs.cram");
        let mut reader = cram::io::Reader::new(Recorder::new(File::open(&path).unwrap()));
        let header = reader.read_header().unwrap();
        let reference = Reference::open(&test_data("contigs.short.fa"), false, 1).unwrap();
        let decoder = CramDecoder::new(Some(reference), &path, &sam::Header::default());
        let mut container = cram::io::reader::Container::default();
        read_container(&mut reader, &mut container).unwrap();
        let error = decoder
            .decode(&container, &header, reader.get_ref().bytes())
            .unwrap_err();
        assert!(error
            .to_string()
            .starts_with("failed to decode a CRAM container of"));
    }

    #[test]
    fn test_quietly_returns_the_value_or_none_on_a_panic() {
        assert_eq!(quietly(|| 7), Some(7));
        assert_eq!(quietly(|| -> u8 { panic!("noodles") }), None);
        assert!(!QUIET.with(std::cell::Cell::get));
    }

    #[test]
    fn test_unchecked_reference_accepts_any_header() {
        let text = "@SQ\tSN:chrM\tLN:1\tM5:x\n";
        let header = sam::io::Reader::new(text.as_bytes()).read_header().unwrap();
        Reference::open(&test_data("ref.fa"), false, 1)
            .unwrap()
            .check_header(&header, Path::new("in.cram"))
            .unwrap();
    }

    #[test]
    fn test_a_missing_fasta_or_fai_is_named() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("ref.fa");
        let err = format!("{:#}", Reference::open(&path, true, 1).err().unwrap());
        assert!(
            err.starts_with("failed to open the --reference FASTA"),
            "{err}"
        );
        assert!(!err.contains(".fai"), "{err}");
        std::fs::write(&path, ">a\nACGT\n").unwrap();
        let err = format!("{:#}", Reference::open(&path, true, 1).err().unwrap());
        assert!(err.starts_with("failed to read the .fai index"), "{err}");
        assert!(err.contains("ref.fa.fai of the --reference FASTA"), "{err}");
        assert!(err.contains("samtools faidx"), "{err}");
    }

    #[test]
    fn test_a_missing_gzi_is_named() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("ref.fa.gz");
        let mut writer = bgzf::io::Writer::new(File::create(&path).unwrap());
        io::Write::write_all(&mut writer, b">a\nACGT\n").unwrap();
        writer.finish().unwrap();
        std::fs::write(dir.path().join("ref.fa.gz.fai"), "a\t4\t3\t4\t5\n").unwrap();
        let err = format!("{:#}", Reference::open(&path, true, 1).err().unwrap());
        assert!(err.starts_with("failed to read the .gzi index"), "{err}");
        assert!(
            err.contains("ref.fa.gz.gzi of the bgzipped --reference FASTA"),
            "{err}"
        );
        assert!(err.contains("bgzip -r"), "{err}");
        bgzf::gzi::fs::write(
            dir.path().join("ref.fa.gz.gzi"),
            &bgzf::gzi::Index::default(),
        )
        .unwrap();
        let reference = Reference::open(&path, true, 1).unwrap();
        assert_eq!(&reference.sequence(b"a").unwrap()[..], b"ACGT");
    }

    #[test]
    fn test_duplicate_contig_names_are_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("dup.fa");
        std::fs::write(&path, ">a\nACGT\n>a\nTTTT\n").unwrap();
        std::fs::write(
            dir.path().join("dup.fa.fai"),
            "a\t4\t3\t4\t5\na\t4\t11\t4\t5\n",
        )
        .unwrap();
        for check in [true, false] {
            let err = Reference::open(&path, check, 1).err().unwrap().to_string();
            assert!(
                err.contains("contig `a` is listed more than once in the .fai index"),
                "{err}"
            );
            assert!(err.contains("dup.fa.fai"), "{err}");
        }
    }

    #[test]
    fn test_reference_reads_a_whole_sequence() {
        let reference = Reference::open(&test_data("contigs.fa"), false, 1).unwrap();
        assert_eq!(reference.sequence(b"ctgB").unwrap().len(), 200);
        assert!(reference.sequence(b"missing").is_err());
    }
}
