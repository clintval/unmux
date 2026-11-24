//! Tag sets: the finite, known set of barcodes a group matches against. A set
//! comes from a `--group NAME=FILE` (a headered, tab-delimited file) or an
//! inline `{...}` declaration, and is the closed set that selector tokens
//! (`group::id-or-seq`) and sheet cells resolve against.
//!
//! A tag file is a headered TSV with `#` comments skipped, an id column (header
//! `id`, `tag`, `name`, `sample_id`, or `sample`) and a sequence column (header
//! `seq`, `sequence`, or `barcode`), plus an optional `sub_sample` column
//! consumed by `--sample-from-group` to set read-group `LB`. Any other column is
//! ignored unless an `--extract NAME=@group::COLUMN` reads it.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::grammar::InlineTag;

/// One tag of a set: an id, its sequence, and an optional sub_sample (`LB`)
/// carried from a file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TagEntry {
    /// The tag id (matches a selector token; doubles as the `SM` for
    /// `--sample-from-group`).
    pub id: String,
    /// The tag sequence (DNA / IUPAC).
    pub seq: String,
    /// An optional sub_sample id from the file's `sub_sample` column.
    pub sub_sample: Option<String>,
}

/// A group's tag set, in declaration / file order.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct TagSet {
    /// The tags, in order.
    pub entries: Vec<TagEntry>,
}

impl TagSet {
    /// Build a tag set from an inline `{...}` declaration (no sub_samples).
    pub fn from_inline(tags: &[InlineTag]) -> Self {
        Self {
            entries: tags
                .iter()
                .map(|t| TagEntry {
                    id: t.id.clone(),
                    seq: t.seq.clone(),
                    sub_sample: None,
                })
                .collect(),
        }
    }

    /// The tag sequences as bytes, in order, for handing to the matcher.
    pub fn seqs(&self) -> Vec<Vec<u8>> {
        self.entries
            .iter()
            .map(|t| t.seq.as_bytes().to_vec())
            .collect()
    }

    /// Resolve a selector token (a tag id or a tag sequence) to its entry
    /// index, or `None` if the token names neither. An id match is preferred
    /// over a sequence match.
    pub fn resolve(&self, token: &str) -> Option<usize> {
        if let Some(i) = self.entries.iter().position(|t| t.id == token) {
            return Some(i);
        }
        self.entries.iter().position(|t| t.seq == token)
    }
}

/// The values of the tag-file columns that `@group::COLUMN` extracts read, keyed
/// by the column as the extract spells it, each holding one value per tag.
pub type TagColumns = HashMap<String, Vec<String>>;

/// Load a group's tag set from a headered, tab-delimited file.
pub fn load_tag_file(path: &Path) -> Result<TagSet> {
    load_tag_file_columns(path, &[]).map(|(tag_set, _)| tag_set)
}

/// Load a group's tag set along with the values of its `wanted` columns.
pub fn load_tag_file_columns(path: &Path, wanted: &[&str]) -> Result<(TagSet, TagColumns)> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read tag file: {}", path.display()))?;
    parse_tag_table(&text, wanted)
        .with_context(|| format!("failed to parse tag file: {}", path.display()))
}

/// Whether a column reference is the tag ID.
pub(crate) fn is_id_column(column: &str) -> bool {
    column.eq_ignore_ascii_case("id")
}

/// Whether a column reference is the tag sequence.
pub(crate) fn is_seq_column(column: &str) -> bool {
    column.eq_ignore_ascii_case("seq")
}

/// Whether `value` can be a SAM tag value with no spaces: every SAM/BAM/CRAM
/// `Z` tag must be printable ASCII.
pub(crate) fn is_column_value(value: &str) -> bool {
    value.bytes().all(|byte| byte.is_ascii_graphic())
}

/// Parse the text of a tag table (headered TSV; blank lines and `#` comments
/// skipped), keeping the values of its `wanted` columns.
fn parse_tag_table(text: &str, wanted: &[&str]) -> Result<(TagSet, TagColumns)> {
    let mut rows = text
        .lines()
        .map(|line| line.trim_end_matches('\r'))
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty() && !line.trim_start().starts_with('#'));

    let (_, header) = rows
        .next()
        .context("tag file is empty (a header row is required)")?;
    let columns: Vec<String> = header
        .split('\t')
        .map(|c| c.trim().to_lowercase())
        .collect();

    // `id` is the generic primary; `sample_id`/`sample` are accepted so a
    // sample-demux sheet (e.g. a `--sample-metadata`-style sheet whose id
    // column is `sample_id`) loads directly as a tag set.
    let id_col = find_column(&columns, &["id", "tag", "name", "sample_id", "sample"])?
        .context("tag file has no id column (expected one of: id, tag, name, sample_id, sample)")?;
    let seq_col = find_column(&columns, &["seq", "sequence", "barcode"])?
        .context("tag file has no sequence column (expected one of: seq, sequence, barcode)")?;
    let sub_sample_col = find_column(&columns, &["sub_sample"])?;
    let mut values = wanted
        .iter()
        .map(|&name| {
            let index = find_wanted_column(header, name, id_col, seq_col)?;
            Ok((name, index, Vec::new()))
        })
        .collect::<Result<Vec<(&str, usize, Vec<String>)>>>()?;

    let mut entries = Vec::new();
    for (line_no, (file_line, row)) in rows.enumerate() {
        let cells: Vec<&str> = row.split('\t').collect();
        let id = cell(&cells, id_col)
            .with_context(|| format!("tag row {} is missing its id column", line_no + 1))?;
        let seq = cell(&cells, seq_col)
            .with_context(|| format!("tag row {} is missing its sequence column", line_no + 1))?;
        if id.is_empty() || seq.is_empty() {
            bail!("tag row {} has an empty id or sequence", line_no + 1);
        }
        let sub_sample = sub_sample_col
            .and_then(|c| cell(&cells, c))
            .filter(|s| !s.is_empty())
            .map(str::to_string);
        for (name, index, column) in &mut values {
            let value = match cell(&cells, *index).filter(|value| !value.is_empty()) {
                Some(value) if !is_column_value(value) => bail!(
                    "line {}: value {value:?} in column `{name}` may hold only printable ASCII \
                     with no whitespace",
                    file_line + 1
                ),
                Some(value) => value,
                None => seq,
            };
            column.push(value.to_string());
        }
        entries.push(TagEntry {
            id: id.to_string(),
            seq: seq.to_string(),
            sub_sample,
        });
    }
    if entries.is_empty() {
        bail!("tag file has a header but no tag rows");
    }
    let values = values
        .into_iter()
        .map(|(name, _, column)| (name.to_string(), column))
        .collect();
    Ok((TagSet { entries }, values))
}

/// Find the column an `@group::COLUMN` extract reads: `id` and `seq` are the id
/// and sequence columns whatever their headers, and any other name matches a
/// header regardless of case.
fn find_wanted_column(header: &str, name: &str, id_col: usize, seq_col: usize) -> Result<usize> {
    if is_id_column(name) {
        return Ok(id_col);
    }
    if is_seq_column(name) {
        return Ok(seq_col);
    }
    let columns: Vec<&str> = header.split('\t').map(str::trim).collect();
    let mut matches = columns
        .iter()
        .enumerate()
        .filter(|(_, column)| column.to_lowercase() == name.to_lowercase());
    let Some((index, _)) = matches.next() else {
        let known: Vec<&str> = columns.into_iter().filter(|c| !c.is_empty()).collect();
        bail!("no column `{name}`; the columns are: {}", known.join(", "));
    };
    if matches.next().is_some() {
        bail!("more than one column `{name}`");
    }
    Ok(index)
}

/// Find the one column whose lowercased header matches an accepted alias.
/// `Ok(None)` when none match (the caller decides if the column is required);
/// an error when two or more columns match aliases for the same field, since an
/// ambiguous header (e.g. both `id` and `tag`) would otherwise silently take
/// whichever came first.
fn find_column(columns: &[String], aliases: &[&str]) -> Result<Option<usize>> {
    let mut matches = columns
        .iter()
        .enumerate()
        .filter(|(_, c)| aliases.contains(&c.as_str()));
    let Some((index, first)) = matches.next() else {
        return Ok(None);
    };
    if let Some((_, second)) = matches.next() {
        bail!(
            "tag file header has more than one column for the same field ({first:?} and {second:?}); \
             use exactly one of: {}",
            aliases.join(", ")
        );
    }
    Ok(Some(index))
}

/// The trimmed value of a cell, or `None` if the row has no such column.
fn cell<'a>(cells: &[&'a str], index: usize) -> Option<&'a str> {
    cells.get(index).map(|c| c.trim())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tag_table_basic_id_seq() {
        let table = "id\tseq\ncbt01\tAACGT\ncbt02\tACAGG\n";
        let set = parse_tag_table(table, &[]).unwrap().0;
        assert_eq!(set.entries.len(), 2);
        assert_eq!(set.entries[0].id, "cbt01");
        assert_eq!(set.entries[0].seq, "AACGT");
        assert_eq!(set.entries[0].sub_sample, None);
    }

    #[test]
    fn test_parse_tag_table_header_aliases() {
        // `name`/`barcode` are accepted aliases for id/seq.
        let table = "name\tbarcode\nt1\tACGT\n";
        let set = parse_tag_table(table, &[]).unwrap().0;
        assert_eq!(set.entries[0].id, "t1");
        assert_eq!(set.entries[0].seq, "ACGT");
    }

    #[test]
    fn test_parse_tag_table_sample_id_alias() {
        // a `--sample-metadata`-style sheet uses `sample_id`/`barcode`; it
        // loads directly as a tag set.
        let table = "sample_id\tbarcode\ns1\tACAGTGGTACTTGATG\n";
        let set = parse_tag_table(table, &[]).unwrap().0;
        assert_eq!(set.entries[0].id, "s1");
        assert_eq!(set.entries[0].seq, "ACAGTGGTACTTGATG");
    }

    #[test]
    fn test_parse_tag_table_skips_comments_and_blanks() {
        let table = "# a comment\n\nid\tseq\n# inner comment\ncbt01\tAACGT\n\n";
        let set = parse_tag_table(table, &[]).unwrap().0;
        assert_eq!(set.entries.len(), 1);
        assert_eq!(set.entries[0].id, "cbt01");
    }

    #[test]
    fn test_parse_tag_table_sub_sample_column() {
        let table = "id\tseq\tsub_sample\ncbt01\tAACGT\tlib01\ncbt02\tACAGG\t\n";
        let set = parse_tag_table(table, &[]).unwrap().0;
        assert_eq!(set.entries[0].sub_sample.as_deref(), Some("lib01"));
        // A blank cell leaves no sub_sample.
        assert_eq!(set.entries[1].sub_sample, None);
    }

    #[test]
    fn test_parse_tag_table_missing_seq_column_is_error() {
        let err = parse_tag_table("id\tnotes\ncbt01\thi\n", &[]).unwrap_err();
        assert!(err.to_string().contains("no sequence column"), "{err}");
    }

    #[test]
    fn test_parse_tag_table_ambiguous_id_header_is_error() {
        // Two columns both alias the id field (`id` and `tag`): rather than
        // silently take the first, the header is rejected so a malformed sheet
        // is caught up front.
        let err = parse_tag_table("id\ttag\tseq\nc1\tc1b\tAACGT\n", &[]).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("more than one column"), "{msg}");
        assert!(
            msg.contains("\"id\"") && msg.contains("\"tag\""),
            "names both: {msg}"
        );
    }

    #[test]
    fn test_parse_tag_table_ambiguous_seq_header_is_error() {
        // `seq` and `barcode` both alias the sequence field.
        let err = parse_tag_table("id\tseq\tbarcode\nc1\tAACGT\tAACGT\n", &[]).unwrap_err();
        assert!(err.to_string().contains("more than one column"), "{err}");
    }

    #[test]
    fn test_parse_tag_table_empty_is_error() {
        let err = parse_tag_table("# only comments\n\n", &[]).unwrap_err();
        assert!(err.to_string().contains("header row is required"), "{err}");
    }

    #[test]
    fn test_parse_tag_table_header_only_is_error() {
        let err = parse_tag_table("id\tseq\n", &[]).unwrap_err();
        assert!(err.to_string().contains("no tag rows"), "{err}");
    }

    #[test]
    fn test_parse_tag_table_crlf_tolerated() {
        let table = "id\tseq\r\ncbt01\tAACGT\r\n";
        let set = parse_tag_table(table, &[]).unwrap().0;
        assert_eq!(set.entries[0].seq, "AACGT");
    }

    #[test]
    fn test_tag_set_resolve_by_id_and_seq() {
        let set = parse_tag_table("id\tseq\ncbt01\tAACGT\ncbt02\tACAGG\n", &[])
            .unwrap()
            .0;
        assert_eq!(set.resolve("cbt01"), Some(0));
        assert_eq!(set.resolve("ACAGG"), Some(1));
        assert_eq!(set.resolve("nope"), None);
    }

    #[test]
    fn test_tag_set_from_inline() {
        let inline = vec![
            InlineTag {
                id: "AACGT".to_string(),
                seq: "AACGT".to_string(),
            },
            InlineTag {
                id: "sci_8".to_string(),
                seq: "ACAGGCG".to_string(),
            },
        ];
        let set = TagSet::from_inline(&inline);
        assert_eq!(set.seqs(), vec![b"AACGT".to_vec(), b"ACAGGCG".to_vec()]);
        assert_eq!(set.resolve("sci_8"), Some(1));
    }

    #[test]
    fn test_load_tag_file_round_trip() {
        use std::io::Write;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        write!(file, "id\tseq\ncbt01\tAACGT\n").unwrap();
        file.flush().unwrap();
        let set = load_tag_file(file.path()).unwrap();
        assert_eq!(set.entries[0].id, "cbt01");
    }

    #[test]
    fn test_parse_tag_table_wanted_columns() {
        let table = "Note\tname\tbarcode\tCell\nx y\tt1\tACGT\tA+B\n\tt2\tTTTT\t \nz\tt3\tGGGG\n";
        let (set, columns) = parse_tag_table(table, &["cell", "id", "NAME", "barcode"]).unwrap();
        assert_eq!(set.entries.len(), 3);
        assert_eq!(columns["cell"], vec!["A+B", "TTTT", "GGGG"]);
        assert_eq!(columns["id"], vec!["t1", "t2", "t3"]);
        assert_eq!(columns["NAME"], vec!["t1", "t2", "t3"]);
        assert_eq!(columns["barcode"], vec!["ACGT", "TTTT", "GGGG"]);
    }

    #[test]
    fn test_parse_tag_table_unknown_column_lists_the_columns() {
        let err = parse_tag_table("id\tseq\tCell\t\nt1\tACGT\tx\n", &["nope"]).unwrap_err();
        assert_eq!(
            err.to_string(),
            "no column `nope`; the columns are: id, seq, Cell"
        );
    }

    #[test]
    fn test_parse_tag_table_repeated_column_is_an_error_only_when_wanted() {
        let table = "id\tseq\tcell\tCell\nt1\tACGT\tx\ty\n";
        assert!(parse_tag_table(table, &[]).is_ok());
        let err = parse_tag_table(table, &["cell"]).unwrap_err();
        assert!(err.to_string().contains("more than one column"), "{err}");
    }

    #[test]
    fn test_parse_tag_table_rejects_unprintable_wanted_values() {
        for bad in ["a b", "a\u{b}b", "a\u{1}", "a\u{7f}", "a\rb", "\u{e9}"] {
            let table = format!("# comment\n\nid\tseq\tcell\nt1\tACGT\tok\nt2\tTTTT\t{bad}\n");
            assert!(parse_tag_table(&table, &[]).is_ok());
            let err = parse_tag_table(&table, &["cell"]).unwrap_err().to_string();
            assert!(err.starts_with("line 5: value"), "{bad:?}: {err}");
        }
    }

    #[test]
    fn test_parse_tag_table_rejects_an_id_with_whitespace_only_when_wanted() {
        let table = "id\tseq\nt 1\tACGT\n";
        assert!(parse_tag_table(table, &[]).is_ok());
        let err = parse_tag_table(table, &["id"]).unwrap_err().to_string();
        assert!(err.starts_with("line 2: value \"t 1\""), "{err}");
    }
}
