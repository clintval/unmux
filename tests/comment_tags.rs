//! Integration: SAM tags embedded in a FASTX read-name comment (the samtools
//! `fastq -T` convention) are lifted into real SAM tags on conversion, so a
//! tag-bearing FASTQ becomes a valid uBAM. This drives the real binary the way
//! a user would: `unmux --out r1.bam < r1.fastq` (here the FASTQ is piped on
//! stdin via `--in 0=-` and the `.bam` extension selects BAM output; a bare
//! `> r1.bam` redirect would not, since stdout mirrors the input format).
//!
//! Lifting is lenient and per-field: a field that is not a well-formed SAM tag
//! (a bare UMI, a CASAVA string, free text) is skipped, while the valid tags in
//! the same comment still lift.

use assert_cmd::Command;
use noodles::bam;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::Data;
use noodles::sam::alignment::RecordBuf;

/// The bytes of a `Z` (string) tag, or `None` if the tag is absent or not a
/// string.
fn z_tag(data: &Data, tag: &[u8; 2]) -> Option<Vec<u8>> {
    match data.get(tag) {
        Some(Value::String(s)) => Some(s.to_vec()),
        _ => None,
    }
}

#[test]
fn fastq_comment_tags_become_real_sam_tags_in_ubam() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.bam");

    // read1's comment mixes a junk field with two well-formed tags; read2 has no
    // comment at all.
    Command::cargo_bin("unmux")
        .unwrap()
        .args(["--in", "0=-", "--out", out.to_str().unwrap()])
        .write_stdin(
            "@read1\tnotatag\tRX:Z:ACGT\tBC:Z:GGGG\nAAAACCCC\n+\nIIIIIIII\n\
             @read2\nTTTTGGGG\n+\nIIIIIIII\n",
        )
        .assert()
        .success();

    let mut reader = bam::io::Reader::new(std::fs::File::open(&out).unwrap());
    let header = reader.read_header().unwrap();
    let records: Vec<RecordBuf> = reader
        .record_bufs(&header)
        .map(|r| r.expect("decode BAM record"))
        .collect();
    assert_eq!(records.len(), 2, "two records round-trip");

    let r1 = &records[0];
    assert_eq!(r1.name().map(|n| n.to_vec()), Some(b"read1".to_vec()));
    // A valid uBAM record is unmapped.
    assert!(r1.flags().is_unmapped(), "the record is unmapped (uBAM)");
    // The two well-formed tags are lifted into real SAM fields.
    assert_eq!(z_tag(r1.data(), b"RX").as_deref(), Some(&b"ACGT"[..]));
    assert_eq!(z_tag(r1.data(), b"BC").as_deref(), Some(&b"GGGG"[..]));
    // The junk field was skipped, not turned into a tag.
    assert!(
        r1.data().get(b"no").is_none(),
        "the `notatag` junk field must not become a tag"
    );

    let r2 = &records[1];
    assert_eq!(r2.name().map(|n| n.to_vec()), Some(b"read2".to_vec()));
    // A read with no comment gains no lifted tags. (unmux assigns every record a
    // fan-out RG tag, so check for the lifted tags specifically rather than an
    // empty tag set.)
    assert!(
        r2.data().get(b"RX").is_none() && r2.data().get(b"BC").is_none(),
        "a read with no comment must not gain RX/BC tags"
    );
}
