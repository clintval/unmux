Use this to split a pooled Illumina run into per-sample unmapped BAMs when the index reads are in their own files, here with a UMI after the i5 index in I2, as in UDiTaS (Giannoukos et al. 2018), a method for measuring gene-editing outcomes.

**Data:** the first {{reads}} reads of the UDiTaS fig2b MiSeq run, with its sample sheet ({{links}}). Its {{samples}} samples share the i5 index N504 and each has its own i7. I1 is the 8 bp i7 index; I2 is the 8 bp i5 index followed by a 10 bp UMI:

{{index_cmd}}

{{index}}

Build a sample sheet and a tag file per index from the run's `sample_info.csv`:

{{sheet}}

{{sheets_cmd}}

{{sheets}}

{{demux}}

Each index group searches only its index read and corrects one substitution, and a read is assigned when its i7 and i5 match one row of the sheet. The UMI is I2 from base 8 on (counting from 0), written to RX with its qualities in QX.

{{counts_cmd}}

{{counts}}

{{view_cmd}}

{{view}}

{{assigned}} of {{reads}} reads were assigned. To check every call, the generator reimplements the rules of `demultiplex()` in UDiTaS's {{helpers}}: each index exact, or one mismatch from its closest barcode, and the pair listed in the sheet. {{masking}} RX and QX equal I2's bases and qualities from base 8 on, on every record.

> [!NOTE]
> UDiTaS writes the UMI with low-quality bases masked to N. unmux writes the bases as read and their qualities to QX, so a later step can weigh them.

**Compared with {{fqtk_version}}.** fqtk reads each read's barcodes from fixed positions given as read structures, and each sample's i7 and i5 joined in one `barcode` column:

{{fqtk_sheet}}

{{fqtk_head}}

{{fqtk}}

fqtk counts mismatches over the joined barcode, where unmux gives each index a budget of its own, so `--max-mismatches 2` admits one mismatch in each index, and `--min-mismatch-delta 1` takes the best sample when it beats the next by one; with `--umi-in-name false` the read names stay as they are, so each read can be matched with unmux's. Per read, against unmux:

{{fqtk_table}}

With these settings fqtk makes unmux's call on all {{fqtk_reads}} reads. They would also admit two mismatches in one index, which unmux does not{{twice}} With fqtk's defaults, one mismatch and a lead of two, the {{lost}} reads with a mismatch in each index are left unmatched.

How unmux expresses fqtk's settings:

| fqtk | unmux |
| --- | --- |
| `--read-structures`, with `B` for the barcode bases of each read | a group's `loc`, such as `i7::loc=2:0:8` |
| `--sample-metadata`, each sample's barcodes joined | `--sample-sheet`, a column per group, and a tag file per group |
| `--max-mismatches`, over all the barcode bases | `dist`, for each group |
| `--min-mismatch-delta` | `mode=nearest`, with `delta` |
| `M` segments, the UMI | `--extract` and `--tag "RX=umi"` |
| `T` segments, written to one FASTQ per sample and read | `--template` and `--out` |
| `--unmatched-prefix` | `--unassigned` |
| `demux-metrics.txt` | `--metrics-per-sample` and `--metrics-summary` |
| IUPAC codes in expected barcodes | IUPAC codes in tags |
| `--header-format illumina` and `--umi-in-name` | no equivalent |

What fqtk does that unmux does not: it can write FASTQ headers in the bcl2fastq style, with the sample barcode in the comment and the UMI in the read name. unmux writes tags such as BC and RX as SAM tags instead, in the record or in the FASTQ comment.
