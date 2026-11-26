//! FASTQ inputs read in chunks of whole records, so that several inputs inflate
//! at once on different workers and any worker parses a chunk.
//!
//! A FASTQ record is four lines, so a chunk is framed by counting line feeds,
//! without parsing. The first input leads: it cuts each chunk at a number of
//! records or of bytes, and every other input then frames the same number of
//! records into its chunk of the same index, so the chunks of one index hold
//! the same fragments. Where an input ends, its chunk keeps the bytes after its
//! last whole record and how it ended, so that parsing it fails or stops
//! exactly where parsing the whole input would.

use std::collections::BTreeMap;
use std::io::{self, BufRead, ErrorKind, Read};
use std::sync::{Mutex, MutexGuard, PoisonError};

/// The chunks the leading input may frame ahead of the slowest other input.
const LEAD: u64 = 4;

/// Whole records framed from one input.
pub(crate) struct Chunk {
    data: Vec<u8>,
    records: usize,
    ending: Option<Ending>,
    /// Where the last whole record lies in `data`.
    last_record: Option<(usize, usize)>,
    /// The input's record before this chunk, if any.
    prior: Vec<u8>,
}

/// How an input ended.
enum Ending {
    /// At the end of its bytes.
    Eof,
    /// With an error from reading or decoding it.
    Failed(io::Error),
}

impl Chunk {
    /// The bytes of the input's record before this chunk, empty before its first.
    pub(crate) fn prior(&self) -> &[u8] {
        &self.prior
    }

    /// A reader of the chunk's bytes that ends as its input did: with the
    /// input's error, if one ended it.
    pub(crate) fn into_reader(self) -> ChunkReader {
        let error = match self.ending {
            Some(Ending::Failed(error)) => Some(error),
            Some(Ending::Eof) | None => None,
        };
        ChunkReader {
            data: self.data,
            at: 0,
            error,
        }
    }
}

/// The bytes of a [`Chunk`], then the error that ended its input, if any.
pub(crate) struct ChunkReader {
    data: Vec<u8>,
    at: usize,
    error: Option<io::Error>,
}

impl BufRead for ChunkReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.at == self.data.len() {
            if let Some(error) = self.error.take() {
                return Err(error);
            }
        }
        Ok(&self.data[self.at..])
    }

    fn consume(&mut self, amount: usize) {
        self.at = (self.at + amount).min(self.data.len());
    }
}

impl Read for ChunkReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let available = self.fill_buf()?;
        let n = available.len().min(buf.len());
        buf[..n].copy_from_slice(&available[..n]);
        self.consume(n);
        Ok(n)
    }
}

/// Frame up to `records` whole records from `reader`, stopping sooner once the
/// chunk holds `max_bytes` at a multiple of `step` records, or where the input
/// ends. The chunk starts with room for `capacity` bytes.
fn frame(
    reader: &mut dyn BufRead,
    records: usize,
    max_bytes: usize,
    step: usize,
    capacity: usize,
) -> Chunk {
    let mut chunk = Chunk {
        data: Vec::with_capacity(capacity),
        records: 0,
        ending: None,
        last_record: None,
        prior: Vec::new(),
    };
    let mut lines = 0;
    let mut record_start = 0;
    while chunk.records < records {
        let buf = match reader.fill_buf() {
            Ok([]) => {
                chunk.ending = Some(Ending::Eof);
                break;
            }
            Ok(buf) => buf,
            Err(error) if error.kind() == ErrorKind::Interrupted => continue,
            Err(error) => {
                chunk.ending = Some(Ending::Failed(error));
                break;
            }
        };
        let mut take = buf.len();
        let mut full = false;
        for at in memchr::memchr_iter(b'\n', buf) {
            lines += 1;
            if lines < 4 {
                continue;
            }
            lines = 0;
            chunk.records += 1;
            let size = chunk.data.len() + at + 1;
            chunk.last_record = Some((record_start, size));
            record_start = size;
            let cut = chunk.records.is_multiple_of(step) && size >= max_bytes;
            if chunk.records == records || cut {
                take = at + 1;
                full = true;
                break;
            }
        }
        chunk.data.extend_from_slice(&buf[..take]);
        reader.consume(take);
        if full {
            break;
        }
    }
    chunk
}

/// The chunks a read completed, each with one chunk per input, by index.
pub(crate) struct Framed {
    pub chunks: Vec<(u64, Vec<Chunk>)>,
    /// The number of chunk indexes, once the leading input has ended.
    pub total: Option<u64>,
}

/// Every FASTQ input of a run, read in chunks.
pub(crate) struct FastqSource {
    inputs: Vec<Mutex<Box<dyn BufRead + Send>>>,
    /// Records to a fragment in each input.
    step: usize,
    /// The most records in a chunk.
    max_records: usize,
    /// The bytes at which the leading input cuts a chunk short.
    max_bytes: usize,
    state: Mutex<Framing>,
}

/// Where each input is in framing its chunks.
struct Framing {
    /// The chunk index each input frames next.
    next: Vec<u64>,
    /// Inputs that have ended, or that must not be read again.
    stopped: Vec<bool>,
    /// The records of each leading chunk that another input has yet to frame,
    /// and whether the leading input ended in it.
    counts: BTreeMap<u64, (usize, bool)>,
    /// Chunks that wait for the other inputs' chunks of their index.
    waiting: BTreeMap<u64, Vec<Option<Chunk>>>,
    /// The index of the leading input's last chunk.
    last: Option<u64>,
    /// The bytes of each input's last chunk, the room to start its next with.
    sizes: Vec<usize>,
    /// The last whole record each input has framed.
    priors: Vec<Vec<u8>>,
}

impl FastqSource {
    /// Read `inputs`, `step` records to a fragment in each, in chunks of up to
    /// `max_records` records that the leading input cuts short at `max_bytes`.
    pub(crate) fn new(
        inputs: Vec<Box<dyn BufRead + Send>>,
        step: usize,
        max_records: usize,
        max_bytes: usize,
    ) -> Self {
        let count = inputs.len();
        Self {
            inputs: inputs.into_iter().map(Mutex::new).collect(),
            step,
            max_records: max_records.max(step) / step * step,
            max_bytes,
            state: Mutex::new(Framing {
                next: vec![0; count],
                stopped: vec![false; count],
                counts: BTreeMap::new(),
                waiting: BTreeMap::new(),
                last: None,
                sizes: vec![0; count],
                priors: vec![Vec::new(); count],
            }),
        }
    }

    fn lock(&self) -> MutexGuard<'_, Framing> {
        self.state.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// The number of inputs.
    pub(crate) fn inputs(&self) -> usize {
        self.inputs.len()
    }

    /// Whether `input` has a chunk to frame now: the leading input while it is
    /// at most [`LEAD`] chunks ahead of every other, and another input once
    /// the leading one has framed its next chunk.
    pub(crate) fn can_read(&self, input: usize) -> bool {
        let state = self.lock();
        let next = state.next[input];
        if state.stopped[input] || state.last.is_some_and(|last| next > last) {
            return false;
        }
        if input > 0 {
            return next < state.next[0];
        }
        let slowest = state.next[1..].iter().min();
        slowest.is_none_or(|&slowest| next < slowest + LEAD)
    }

    /// Stop reading every input.
    pub(crate) fn stop(&self) {
        self.lock().stopped.fill(true);
    }

    /// Frame the next chunk of `input`, returning the chunk indexes it
    /// completes. Only one worker reads an input at a time.
    pub(crate) fn read(&self, input: usize) -> Framed {
        let (index, records, max_bytes, step, capacity, prior) = {
            let mut state = self.lock();
            let index = state.next[input];
            let capacity = state.sizes[input] + state.sizes[input] / 8;
            let prior = std::mem::take(&mut state.priors[input]);
            if input == 0 {
                let (records, bytes) = (self.max_records, self.max_bytes);
                (index, records, bytes, self.step, capacity, prior)
            } else {
                // Past the leading input's last whole record, frame two more
                // records: its bytes after them may hold one more record.
                let (count, ended) = state.counts.get(&index).copied().unwrap_or((0, true));
                let records = if ended { count + 2 } else { count };
                (index, records, usize::MAX, 1, capacity, prior)
            }
        };
        let mut chunk = {
            let mut reader = self.inputs[input]
                .lock()
                .unwrap_or_else(PoisonError::into_inner);
            frame(reader.as_mut(), records, max_bytes, step, capacity)
        };
        let mut state = self.lock();
        let ended = chunk.ending.is_some();
        state.next[input] = index + 1;
        state.stopped[input] |= ended;
        state.sizes[input] = chunk.data.len();
        state.priors[input] = match chunk.last_record {
            Some((start, end)) => chunk.data[start..end].to_vec(),
            None => prior.clone(),
        };
        chunk.prior = prior;
        let mut framed = Framed {
            chunks: Vec::new(),
            total: None,
        };
        if input == 0 {
            if self.inputs.len() > 1 {
                state.counts.insert(index, (chunk.records, ended));
            }
            if ended {
                state.last = Some(index);
                framed.total = Some(index + 1);
            }
        }
        let inputs = self.inputs.len();
        let slots = state
            .waiting
            .entry(index)
            .or_insert_with(|| (0..inputs).map(|_| None).collect());
        slots[input] = Some(chunk);
        if slots.iter().all(Option::is_some) {
            let slots = state.waiting.remove(&index).unwrap_or_default();
            framed
                .chunks
                .push((index, slots.into_iter().flatten().collect()));
        }
        if let Some(&slowest) = state.next[1..].iter().min() {
            state.counts.retain(|&at, _| at >= slowest);
        }
        framed
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fastq(names: impl IntoIterator<Item = usize>) -> Vec<u8> {
        names
            .into_iter()
            .flat_map(|i| format!("@r{i}\nACGT\n+\nIIII\n").into_bytes())
            .collect()
    }

    fn source(
        inputs: Vec<Vec<u8>>,
        step: usize,
        max_records: usize,
        max_bytes: usize,
    ) -> FastqSource {
        let inputs = inputs
            .into_iter()
            .map(|bytes| Box::new(io::Cursor::new(bytes)) as Box<dyn BufRead + Send>)
            .collect();
        FastqSource::new(inputs, step, max_records, max_bytes)
    }

    /// Read every input whenever it can be, as the pool would, collecting the
    /// completed chunk indexes in order with the total.
    fn drain(source: &FastqSource) -> (BTreeMap<u64, Vec<Chunk>>, Option<u64>) {
        let mut chunks = BTreeMap::new();
        let mut total = None;
        loop {
            let readable = (0..source.inputs()).find(|&input| source.can_read(input));
            let Some(input) = readable else {
                return (chunks, total);
            };
            let framed = source.read(input);
            chunks.extend(framed.chunks);
            total = total.or(framed.total);
        }
    }

    fn bytes(chunk: Chunk) -> Vec<u8> {
        let mut out = Vec::new();
        chunk.into_reader().read_to_end(&mut out).unwrap();
        out
    }

    #[test]
    fn test_chunks_hold_whole_records_and_cover_the_input() {
        let input = fastq(0..10);
        let (chunks, total) = drain(&source(vec![input.clone()], 1, 3, usize::MAX));
        assert_eq!(total, Some(4));
        let records: Vec<usize> = chunks.values().map(|c| c[0].records).collect();
        assert_eq!(records, [3, 3, 3, 1]);
        let joined: Vec<u8> = chunks
            .into_values()
            .flat_map(|mut c| bytes(c.remove(0)))
            .collect();
        assert_eq!(joined, input);
    }

    #[test]
    fn test_each_chunk_carries_the_record_before_it() {
        let (chunks, _) = drain(&source(vec![fastq(0..7), fastq(0..7)], 1, 3, usize::MAX));
        for input in 0..2 {
            let priors: Vec<&[u8]> = chunks.values().map(|c| c[input].prior()).collect();
            assert_eq!(priors, [&b""[..], &fastq([2]), &fastq([5])]);
        }
    }

    #[test]
    fn test_the_leading_input_cuts_chunks_at_its_byte_limit_on_a_fragment() {
        let record = fastq([0]).len();
        let (chunks, _) = drain(&source(vec![fastq(0..9)], 2, 100, record + 1));
        let records: Vec<usize> = chunks.values().map(|c| c[0].records).collect();
        assert_eq!(records, [2, 2, 2, 2, 1]);
    }

    #[test]
    fn test_every_input_frames_the_leading_input_record_counts() {
        let (chunks, total) = drain(&source(vec![fastq(0..7), fastq(0..7)], 1, 100, 40));
        assert_eq!(total, Some(3));
        let records: Vec<usize> = chunks.values().map(|c| c[1].records).collect();
        assert_eq!(records, [3, 3, 1]);
        for chunk in chunks.values() {
            assert_eq!(chunk[0].data, chunk[1].data);
        }
    }

    #[test]
    fn test_a_longer_other_input_frames_its_records_past_the_end() {
        let (mut chunks, total) = drain(&source(vec![fastq(0..4), fastq(0..9)], 1, 3, usize::MAX));
        assert_eq!(total, Some(2));
        let last = chunks.remove(&1).unwrap();
        assert_eq!(last[0].records, 1);
        assert_eq!(last[1].records, 3);
    }

    #[test]
    fn test_a_shorter_other_input_ends_early() {
        let (chunks, _) = drain(&source(vec![fastq(0..9), fastq(0..4)], 1, 3, usize::MAX));
        let ended: Vec<(u64, usize)> = chunks
            .iter()
            .filter(|(_, c)| c[1].ending.is_some())
            .map(|(&index, c)| (index, c[1].records))
            .collect();
        assert_eq!(ended, [(1, 1)]);
        assert!(chunks.keys().all(|&index| index <= 1));
    }

    #[test]
    fn test_a_partial_last_record_stays_in_the_last_chunk() {
        let mut input = fastq(0..2);
        input.extend_from_slice(b"@r2\nAC");
        let (chunks, total) = drain(&source(vec![input.clone()], 1, 5, usize::MAX));
        assert_eq!(total, Some(1));
        assert_eq!(chunks[&0][0].records, 2);
        let only = chunks.into_values().next().unwrap().remove(0);
        assert_eq!(bytes(only), input);
    }

    /// A reader that yields its bytes and then fails.
    struct Failing(io::Cursor<Vec<u8>>);

    impl Read for Failing {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            match self.0.read(buf)? {
                0 => Err(io::Error::new(
                    ErrorKind::InvalidData,
                    "corrupt deflate stream",
                )),
                n => Ok(n),
            }
        }
    }

    #[test]
    fn test_a_read_error_ends_the_chunk_after_its_bytes() {
        let mut input = fastq(0..2);
        input.extend_from_slice(b"@r2\n");
        let failing = io::BufReader::new(Failing(io::Cursor::new(input.clone())));
        let source = FastqSource::new(vec![Box::new(failing)], 1, 10, usize::MAX);
        let (mut chunks, total) = drain(&source);
        assert_eq!(total, Some(1));
        let mut reader = chunks.remove(&0).unwrap().remove(0).into_reader();
        let mut out = vec![0; input.len()];
        reader.read_exact(&mut out).unwrap();
        assert_eq!(out, input);
        let error = reader.fill_buf().unwrap_err();
        assert_eq!(error.to_string(), "corrupt deflate stream");
    }

    #[test]
    fn test_an_empty_input_is_one_empty_chunk() {
        let (chunks, total) = drain(&source(vec![Vec::new(), Vec::new()], 1, 3, usize::MAX));
        assert_eq!(total, Some(1));
        assert_eq!(chunks[&0].len(), 2);
        assert!(chunks[&0].iter().all(|chunk| chunk.data.is_empty()));
    }
}
