//! `unmux`: flexible record parsing and demultiplexing to FASTX/SAM/BAM/CRAM
//! (SAM/BAM/CRAM written unmapped), splitcode-style, on the `sassy`
//! approximate matcher.
//!
//! The command language is one quote-free grammar shared by the CLI and sheets.
//! The demux command exposes a `DemuxArgs` struct and a `run_demux` entry
//! point, re-exported here for the binary to dispatch.
#![warn(missing_docs)]

pub mod demux;
pub mod extract;
pub mod fanout;
pub mod grammar;
pub mod input;
pub mod iupac;
pub mod matcher;
pub mod metrics;
pub mod output;
pub mod qc;
pub mod tags;
pub mod writer;

pub use demux::{run_demux, DemuxArgs};
pub use extract::{extract, Extracted, MatchSpan, Segment};
pub use fanout::{
    compile_routing, expand_from_group, load_sample_sheet, parse_sample_sheet, Disposition,
    RemoveTarget, Routing, Target,
};
pub use grammar::{parse_demux, DemuxPlan};
pub use input::{sniff_bytes, sniff_input, Fragment, FragmentReader, InputRecord, SniffedFormat};
pub use matcher::{match_group, next_window, CompiledGroup, GroupOutcome, MatchStrand, TagMatch};
pub use tags::{load_tag_file, TagEntry, TagSet};
pub use writer::{output_format, OutputFormat, OutputWriter};
