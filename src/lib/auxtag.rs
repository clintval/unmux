//! The `@tag::XX` split: read a per-record aux-tag value, detect grouping, and
//! guard the high-cardinality output (one file open at a time when grouped, a
//! descriptor-safe all-open fallback otherwise).

use std::collections::HashSet;
use std::hash::{Hash, Hasher};

use anyhow::{bail, Result};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;

use crate::input::Fragment;

/// The value of aux tag `tag` for a fragment: the first segment that carries a
/// string-typed `XX`, or `None` (absent on every segment, or a non-string
/// type, which is treated as absent, matching `samtools split`).
pub fn match_aux_tag(tag: [u8; 2], fragment: &Fragment) -> Option<Vec<u8>> {
    let key = Tag::new(tag[0], tag[1]);
    for record in &fragment.records {
        if let Some(data) = &record.tags {
            if let Some(Value::String(value)) = data.get(&key) {
                return Some(AsRef::<[u8]>::as_ref(value).to_vec());
            }
        }
    }
    None
}

/// Detects a split-key value that reappears after its output file was closed,
/// which means a supposedly grouped input is not actually grouped and a naive
/// reopen would truncate the earlier output. Holds a fingerprint per closed
/// value, so it is safe for high-cardinality streams.
#[derive(Default)]
pub struct ClosedKeys {
    closed: HashSet<u64>,
}

impl ClosedKeys {
    /// Error if `key` was already closed (a recurrence, so the input is not
    /// grouped by the split tag); otherwise a no-op.
    pub fn enter(&mut self, key: &[u8]) -> Result<()> {
        if self.closed.contains(&fingerprint(key)) {
            bail!(
                "value `{}` reappeared after its group was closed, so the input is not grouped by the split tag; sort by the tag first (e.g. `samtools sort -t <TAG>`) and re-run",
                String::from_utf8_lossy(key)
            );
        }
        Ok(())
    }

    /// Record `key` as closed (called when its output file is finalized).
    pub fn close(&mut self, key: &[u8]) {
        self.closed.insert(fingerprint(key));
    }
}

fn fingerprint(key: &[u8]) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    key.hash(&mut hasher);
    hasher.finish()
}

/// Raise this process's soft open-file limit toward the hard cap (a process may
/// do this without privileges) and return the effective soft limit. Best-effort:
/// on any error the current soft limit is returned unchanged.
pub fn raise_open_file_limit() -> u64 {
    let current = rlimit::Resource::NOFILE
        .get()
        .map(|(soft, _hard)| soft)
        .unwrap_or(0);
    rlimit::increase_nofile_limit(u64::MAX).unwrap_or(current)
}

/// How many output files may be open at once: the soft limit minus a reserve
/// for stdio, the input files, and slack. Floors at zero.
pub fn open_file_ceiling(soft_limit: u64, reserved: u64) -> usize {
    soft_limit.saturating_sub(reserved) as usize
}

/// Whether an I/O error is the OS refusing another open file
/// (`EMFILE`, per-process, or `ENFILE`, system-wide).
pub fn is_too_many_open_files(err: &std::io::Error) -> bool {
    matches!(err.raw_os_error(), Some(libc::EMFILE) | Some(libc::ENFILE))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::{Fragment, InputRecord};
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    use noodles::sam::alignment::record_buf::Data;

    fn frag_with(tag: &[u8; 2], value: Value) -> Fragment {
        let mut data = Data::default();
        data.insert(Tag::new(tag[0], tag[1]), value);
        Fragment {
            records: vec![InputRecord {
                name: b"r1".to_vec(),
                bases: b"ACGT".to_vec(),
                quals: None,
                tags: Some(data),
            }],
        }
    }

    #[test]
    fn reads_a_string_value() {
        let f = frag_with(b"CB", Value::String("AAACCC".into()));
        assert_eq!(match_aux_tag(*b"CB", &f), Some(b"AAACCC".to_vec()));
    }

    #[test]
    fn absent_tag_is_none() {
        let f = frag_with(b"CB", Value::String("AAACCC".into()));
        assert_eq!(match_aux_tag(*b"RX", &f), None);
    }

    #[test]
    fn non_string_value_is_none() {
        let f = frag_with(b"CB", Value::Int32(7));
        assert_eq!(match_aux_tag(*b"CB", &f), None);
    }

    #[test]
    fn entering_a_fresh_key_is_ok_then_closing_marks_it() {
        let mut closed = ClosedKeys::default();
        assert!(closed.enter(b"AAAA").is_ok());
        closed.close(b"AAAA");
        // A different key is still fine.
        assert!(closed.enter(b"CCCC").is_ok());
    }

    #[test]
    fn re_entering_a_closed_key_errors() {
        let mut closed = ClosedKeys::default();
        closed.enter(b"AAAA").unwrap();
        closed.close(b"AAAA");
        let err = closed.enter(b"AAAA").err().unwrap();
        assert!(err.to_string().contains("reappeared"));
    }

    #[test]
    fn ceiling_subtracts_the_reserve_and_floors_at_zero() {
        assert_eq!(open_file_ceiling(1024, 16), 1008);
        assert_eq!(open_file_ceiling(10, 16), 0); // never negative
    }

    #[test]
    fn emfile_and_enfile_are_recognized() {
        let emfile = std::io::Error::from_raw_os_error(libc::EMFILE);
        let enfile = std::io::Error::from_raw_os_error(libc::ENFILE);
        let other = std::io::Error::from_raw_os_error(libc::ENOENT);
        assert!(is_too_many_open_files(&emfile));
        assert!(is_too_many_open_files(&enfile));
        assert!(!is_too_many_open_files(&other));
    }

    #[test]
    fn raising_the_open_file_limit_returns_a_nonzero_soft_limit() {
        assert!(raise_open_file_limit() > 0);
    }
}
