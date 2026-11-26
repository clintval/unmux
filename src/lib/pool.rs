//! One pool of worker threads that runs every stage of a run: reading input
//! batches, processing them, merging their results in input order, and
//! compressing the output blocks those results fill. Any worker takes whichever
//! stage has work, so the threads go where the time goes instead of to fixed
//! roles. Reading a stream and merging are serial; processing and compressing
//! run on every worker with nothing more urgent to do.
//!
//! The queues between the stages are bounded: a stream is read only while few
//! batches are in flight, and a result is merged only while few blocks wait to
//! be compressed, so memory stays bounded whatever the speed of each stage.

use std::collections::{BTreeMap, VecDeque};
use std::sync::{Condvar, Mutex, MutexGuard, PoisonError};

use anyhow::{anyhow, Context, Error, Result};

/// The stages a [`run`] drives. The pool reads each stream on one worker at a
/// time and merges and finishes on one worker at a time; the other stages run
/// on any number of workers at once.
pub(crate) trait Stages: Sync {
    /// A unit of input, read from the streams.
    type Batch: Send;
    /// The result of processing a batch.
    type Done: Send;
    /// A unit of output to compress, filled by merging results.
    type Block: Send;
    /// The state each worker keeps across the jobs it runs.
    type Worker;

    /// The number of input streams.
    fn streams(&self) -> usize;
    /// A new worker's state.
    fn worker(&self) -> Self::Worker;
    /// Whether reading `stream` now would make progress.
    fn can_read(&self, stream: usize) -> bool;
    /// Read from `stream`.
    fn read(&self, worker: &mut Self::Worker, stream: usize) -> Readout<Self::Batch>;
    /// Process a batch.
    fn process(&self, worker: &mut Self::Worker, batch: Self::Batch) -> Self::Done;
    /// Merge the result of the next batch in input order, adding the blocks it
    /// fills to `blocks`.
    fn merge(&self, done: Self::Done, blocks: &mut Vec<Self::Block>) -> Result<()>;
    /// Add the last blocks to `blocks` once every batch is merged.
    fn finish(&self, blocks: &mut Vec<Self::Block>) -> Result<()>;
    /// Compress and write a block.
    fn compress(&self, worker: &mut Self::Worker, block: Self::Block) -> Result<()>;
}

/// What one read produced.
pub(crate) struct Readout<B> {
    /// Batches completed by this read, each with its index in input order.
    pub batches: Vec<(u64, B)>,
    /// The end of the whole input, if this read reached it.
    pub end: Option<End>,
}

impl<B> Readout<B> {
    /// A read that produced nothing.
    pub fn nothing() -> Self {
        Self {
            batches: Vec::new(),
            end: None,
        }
    }
}

/// How the input ended.
pub(crate) enum End {
    /// After this many batches.
    Exhausted(u64),
    /// With an error at this batch index, reported once every batch before it
    /// is merged.
    Failed(u64, Error),
}

/// The most batches read but not yet merged, per worker.
const BATCHES_PER_WORKER: usize = 2;

/// The most blocks waiting to be compressed, per worker, before merging waits.
const BLOCKS_PER_WORKER: usize = 8;

/// A job a worker takes from the pool.
enum Job<S: Stages> {
    Read(usize),
    Process(u64, S::Batch),
    Merge(S::Done),
    Finish,
    Compress(S::Block),
}

/// What the workers share, behind the pool's lock.
struct State<S: Stages> {
    reading: Vec<bool>,
    batches: BTreeMap<u64, S::Batch>,
    processing: usize,
    done: BTreeMap<u64, S::Done>,
    next_merge: u64,
    merging: bool,
    total: Option<u64>,
    read_error: Option<Error>,
    finished: bool,
    blocks: VecDeque<S::Block>,
    compressing: usize,
    error: Option<Error>,
    waiting: usize,
}

impl<S: Stages> State<S> {
    /// Batches read and not yet merged.
    fn in_flight(&self) -> usize {
        self.batches.len() + self.processing + self.done.len() + usize::from(self.merging)
    }

    /// Whether every batch is merged and every block written.
    fn complete(&self) -> bool {
        self.finished && self.blocks.is_empty() && self.compressing == 0
    }

    /// Whether any job is running.
    fn busy(&self) -> bool {
        self.reading.contains(&true) || self.processing > 0 || self.merging || self.compressing > 0
    }

    /// Fail with the read error once every batch before it is merged.
    fn reach_read_error(&mut self) {
        if self.read_error.is_some() && self.total == Some(self.next_merge) && !self.merging {
            let error = self.read_error.take().expect("a read error is pending");
            self.error.get_or_insert(error);
        }
    }
}

/// A pool of workers and the state they share.
struct Pool<'a, S: Stages> {
    stages: &'a S,
    state: Mutex<State<S>>,
    wake: Condvar,
    batch_cap: usize,
    block_cap: usize,
}

impl<S: Stages> Pool<'_, S> {
    fn lock(&self) -> MutexGuard<'_, State<S>> {
        self.state.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// The next job by priority, waiting for one, or `None` once the run is
    /// over or has failed.
    fn next_job(&self) -> Option<Job<S>> {
        let mut state = self.lock();
        loop {
            if state.error.is_some() || state.complete() {
                return None;
            }
            if let Some(job) = self.pick(&mut state) {
                return Some(job);
            }
            if !state.busy() {
                state.error = Some(anyhow!(
                    "internal error: the worker pool stalled before the input was merged"
                ));
                self.wake.notify_all();
                return None;
            }
            state.waiting += 1;
            state = self
                .wake
                .wait(state)
                .unwrap_or_else(PoisonError::into_inner);
            state.waiting -= 1;
        }
    }

    /// Take the most urgent job. Merging and reading are serial, so they come
    /// first to keep the stages behind them busy; compressing comes before
    /// processing so that output drains before more is made.
    fn pick(&self, state: &mut State<S>) -> Option<Job<S>> {
        if !state.merging && state.blocks.len() < self.block_cap {
            if let Some(done) = state.done.remove(&state.next_merge) {
                state.merging = true;
                return Some(Job::Merge(done));
            }
            let merged_all = state.total == Some(state.next_merge);
            if merged_all && state.read_error.is_none() && !state.finished {
                state.merging = true;
                return Some(Job::Finish);
            }
        }
        if state.in_flight() < self.batch_cap {
            let free = (0..state.reading.len())
                .find(|&stream| !state.reading[stream] && self.stages.can_read(stream));
            if let Some(stream) = free {
                state.reading[stream] = true;
                return Some(Job::Read(stream));
            }
        }
        if let Some(block) = state.blocks.pop_front() {
            state.compressing += 1;
            return Some(Job::Compress(block));
        }
        if let Some((index, batch)) = state.batches.pop_first() {
            state.processing += 1;
            return Some(Job::Process(index, batch));
        }
        None
    }

    /// Run one job and publish what it produced.
    fn run(&self, worker: &mut S::Worker, job: Job<S>) {
        match job {
            Job::Read(stream) => {
                let readout = self.stages.read(worker, stream);
                let mut state = self.lock();
                state.reading[stream] = false;
                match readout.end {
                    Some(End::Exhausted(total)) => state.total = Some(total),
                    Some(End::Failed(index, error))
                        if state.total.is_none_or(|total| index < total) =>
                    {
                        state.total = Some(index);
                        state.batches.retain(|&at, _| at < index);
                        state.read_error = Some(error);
                        state.reach_read_error();
                    }
                    Some(End::Failed(..)) | None => {}
                }
                let total = state.total;
                let batches = readout.batches.into_iter();
                let kept = batches.filter(|&(at, _)| total.is_none_or(|total| at < total));
                state.batches.extend(kept);
                self.notify(&state);
            }
            Job::Process(index, batch) => {
                let done = self.stages.process(worker, batch);
                let mut state = self.lock();
                state.processing -= 1;
                state.done.insert(index, done);
                self.notify(&state);
            }
            Job::Merge(done) => {
                let mut blocks = Vec::new();
                let merged = self.stages.merge(done, &mut blocks);
                let mut state = self.lock();
                state.merging = false;
                state.blocks.extend(blocks);
                match merged {
                    Ok(()) => {
                        state.next_merge += 1;
                        state.reach_read_error();
                    }
                    Err(error) => {
                        state.error.get_or_insert(error);
                    }
                }
                self.notify(&state);
            }
            Job::Finish => {
                let mut blocks = Vec::new();
                let finished = self.stages.finish(&mut blocks);
                let mut state = self.lock();
                state.merging = false;
                state.blocks.extend(blocks);
                match finished {
                    Ok(()) => state.finished = true,
                    Err(error) => {
                        state.error.get_or_insert(error);
                    }
                }
                self.notify(&state);
            }
            Job::Compress(block) => {
                let written = self.stages.compress(worker, block);
                let mut state = self.lock();
                state.compressing -= 1;
                if let Err(error) = written {
                    state.error.get_or_insert(error);
                }
                self.notify(&state);
            }
        }
    }

    /// Wake the waiting workers to look for new work.
    fn notify(&self, state: &State<S>) {
        if state.waiting > 0 {
            self.wake.notify_all();
        }
    }

    /// Fail the run and wake every worker, as a worker unwinds from a panic.
    fn abandon(&self) {
        let mut state = self.lock();
        state
            .error
            .get_or_insert_with(|| anyhow!("internal error: a worker thread panicked"));
        self.wake.notify_all();
    }
}

/// Fails the run if the worker that holds it unwinds.
struct Unwind<'a, 'b, S: Stages>(&'a Pool<'b, S>);

impl<S: Stages> Drop for Unwind<'_, '_, S> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            self.0.abandon();
        }
    }
}

/// The stack of each worker; decoding CRAM 3.1 recurses deeper than the
/// default 2 MiB allows in optimized builds.
const WORKER_STACK_BYTES: usize = 64 << 20;

/// Run every stage of `stages` on `workers` threads until the input is merged
/// and every block written, or until the first error in input order.
pub(crate) fn run<S: Stages>(stages: &S, workers: usize) -> Result<()> {
    let workers = workers.max(1);
    let pool = Pool {
        stages,
        state: Mutex::new(State {
            reading: vec![false; stages.streams()],
            batches: BTreeMap::new(),
            processing: 0,
            done: BTreeMap::new(),
            next_merge: 0,
            merging: false,
            total: None,
            read_error: None,
            finished: false,
            blocks: VecDeque::new(),
            compressing: 0,
            error: None,
            waiting: 0,
        }),
        wake: Condvar::new(),
        batch_cap: BATCHES_PER_WORKER * workers + 2,
        block_cap: BLOCKS_PER_WORKER * workers,
    };
    std::thread::scope(|scope| -> Result<()> {
        for index in 0..workers {
            let pool = &pool;
            std::thread::Builder::new()
                .name(format!("unmux-worker-{index}"))
                .stack_size(WORKER_STACK_BYTES)
                .spawn_scoped(scope, move || {
                    let _unwind = Unwind(pool);
                    let mut worker = pool.stages.worker();
                    while let Some(job) = pool.next_job() {
                        pool.run(&mut worker, job);
                    }
                })
                .context("failed to start a worker thread")?;
        }
        Ok(())
    })?;
    let mut state = pool.lock();
    match state.error.take() {
        Some(error) => Err(error),
        None => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::sync::atomic::{AtomicU64, Ordering};

    /// Streams of numbers read in batches of `size`, summed in order into
    /// blocks, with errors at chosen points.
    struct Numbers {
        streams: usize,
        batches: u64,
        size: u64,
        read_error_at: Option<u64>,
        process_error_at: Option<u64>,
        next: Mutex<Vec<u64>>,
        merged: Mutex<Vec<u64>>,
        written: Mutex<Vec<u64>>,
        finished: AtomicU64,
    }

    impl Numbers {
        fn new(streams: usize, batches: u64) -> Self {
            Self {
                streams,
                batches,
                size: 3,
                read_error_at: None,
                process_error_at: None,
                next: Mutex::new(vec![0; streams]),
                merged: Mutex::new(Vec::new()),
                written: Mutex::new(Vec::new()),
                finished: AtomicU64::new(0),
            }
        }
    }

    impl Stages for Numbers {
        type Batch = Vec<u64>;
        type Done = Result<u64>;
        type Block = u64;
        type Worker = u64;

        fn streams(&self) -> usize {
            self.streams
        }

        fn worker(&self) -> u64 {
            0
        }

        fn can_read(&self, stream: usize) -> bool {
            let next = self.next.lock().unwrap();
            let slowest = next.iter().min().copied().unwrap();
            let limit = match self.read_error_at {
                Some(at) => at + u64::from(stream == 0),
                None => u64::MAX,
            };
            next[stream] < self.batches.max(1).min(limit) && next[stream] <= slowest + 1
        }

        fn read(&self, worker: &mut u64, stream: usize) -> Readout<Vec<u64>> {
            *worker += 1;
            let mut next = self.next.lock().unwrap();
            let index = next[stream];
            next[stream] += 1;
            let mut readout = Readout::nothing();
            if self.read_error_at == Some(index) {
                readout.end = Some(End::Failed(index, anyhow!("bad read at {index}")));
            } else if self.batches == 0 {
                readout.end = Some(End::Exhausted(0));
            } else if next.iter().all(|&n| n > index) {
                let first = index * self.size;
                readout
                    .batches
                    .push((index, (first..first + self.size).collect()));
                if index + 1 == self.batches {
                    readout.end = Some(End::Exhausted(self.batches));
                }
            }
            readout
        }

        fn process(&self, _: &mut u64, batch: Vec<u64>) -> Result<u64> {
            let index = batch[0] / self.size;
            if self.process_error_at == Some(index) {
                return Err(anyhow!("bad record in {index}"));
            }
            Ok(batch.iter().sum())
        }

        fn merge(&self, done: Result<u64>, blocks: &mut Vec<u64>) -> Result<()> {
            let sum = done?;
            self.merged.lock().unwrap().push(sum);
            blocks.push(sum);
            Ok(())
        }

        fn finish(&self, blocks: &mut Vec<u64>) -> Result<()> {
            self.finished.fetch_add(1, Ordering::Relaxed);
            blocks.push(u64::MAX);
            Ok(())
        }

        fn compress(&self, _: &mut u64, block: u64) -> Result<()> {
            self.written.lock().unwrap().push(block);
            Ok(())
        }
    }

    /// The sums of each batch of `numbers`, in order.
    fn sums(numbers: &Numbers, batches: u64) -> Vec<u64> {
        (0..batches)
            .map(|index| (0..numbers.size).map(|i| index * numbers.size + i).sum())
            .collect()
    }

    #[test]
    fn test_results_merge_in_input_order_at_every_worker_count() {
        for workers in [1, 2, 3, 8] {
            for streams in [1, 2, 4] {
                let numbers = Numbers::new(streams, 200);
                run(&numbers, workers).unwrap();
                assert_eq!(*numbers.merged.lock().unwrap(), sums(&numbers, 200));
                let mut written = numbers.written.lock().unwrap().clone();
                written.sort_unstable();
                let mut expected = sums(&numbers, 200);
                expected.push(u64::MAX);
                expected.sort_unstable();
                assert_eq!(written, expected);
                assert_eq!(numbers.finished.load(Ordering::Relaxed), 1);
            }
        }
    }

    #[test]
    fn test_an_empty_input_still_finishes() {
        for streams in [1, 3] {
            let empty = Numbers::new(streams, 0);
            run(&empty, 4).unwrap();
            assert!(empty.merged.lock().unwrap().is_empty());
            assert_eq!(*empty.written.lock().unwrap(), [u64::MAX]);
            assert_eq!(empty.finished.load(Ordering::Relaxed), 1);
        }
    }

    #[test]
    fn test_a_read_error_follows_every_batch_before_it() {
        for workers in [1, 2, 8] {
            let mut numbers = Numbers::new(2, 50);
            numbers.read_error_at = Some(20);
            let error = run(&numbers, workers).unwrap_err();
            assert_eq!(error.to_string(), "bad read at 20");
            assert_eq!(*numbers.merged.lock().unwrap(), sums(&numbers, 20));
            assert_eq!(numbers.finished.load(Ordering::Relaxed), 0);
        }
    }

    #[test]
    fn test_the_first_error_in_input_order_wins() {
        for workers in [1, 2, 8] {
            let mut numbers = Numbers::new(1, 50);
            numbers.process_error_at = Some(10);
            numbers.read_error_at = Some(30);
            let error = run(&numbers, workers).unwrap_err();
            assert_eq!(error.to_string(), "bad record in 10");
            assert_eq!(*numbers.merged.lock().unwrap(), sums(&numbers, 10));
        }
    }

    #[test]
    fn test_a_read_error_before_any_batch_fails_the_run() {
        let mut numbers = Numbers::new(1, 50);
        numbers.read_error_at = Some(0);
        let error = run(&numbers, 3).unwrap_err();
        assert_eq!(error.to_string(), "bad read at 0");
        assert!(numbers.merged.lock().unwrap().is_empty());
    }
}
