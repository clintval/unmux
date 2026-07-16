//! The `@tag::XX` split: read a per-record aux-tag value, detect grouping, and
//! guard the high-cardinality output (one file open at a time when grouped, a
//! descriptor-safe all-open fallback otherwise).

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
}
