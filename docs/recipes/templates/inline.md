Use this when the sample barcode is the first bases of the read rather than an index read, as in RAD-seq: unmux checks the restriction site after it, demultiplexes on the barcode, trims it off, and keeps it in BC with its qualities in QT.

**Data:** 250 reads of each of the 32 RAD-seq runs of ground beetles in study {{study_link}} (Takahashi et al. 2014), one run per sample, fetched with fg-sra and pooled one after another, so the sample of every read is known. Each 101 bp read starts with the sample's 5 bp barcode and the PstI site remnant TGCAG; the barcodes are the ones each run's BioSample record (SAMD00008878 to SAMD00008909) gives, at least 3 substitutions apart.

{{fetch}}

{{head_cmd}}

{{head}}

{{demux}}

The `site` group requires the remnant right after the barcode, so a read without it would be unassigned; the insert keeps the remnant, as RAD pipelines expect.

A read whose barcode had one substitution:

{{view_cmd}}

{{view}}

All {{reads}} reads went to their run's sample and none to a wrong one; {{site}} carry the remnant exactly. On every record BC and QT equal the read's first 5 bases and qualities, and the record holds the rest of the read.

> [!NOTE]
> The submitter deposited these runs already split by barcode, allowing one substitution, so no read here has more errors than `dist=1` corrects.

> [!TIP]
> With `dist=0`, {{exact}} reads are assigned and the {{inexact}} with a substitution in the barcode are not.

**Compared with {{fqtk_version}}.** The same pool through fqtk, with one mismatch, as `dist=1` allows, and a lead of one over the next barcode; {{ref:uditas}} gives how its settings map to unmux's:

{{fqtk}}

{{fqtk_table}}

fqtk makes unmux's call on {{fqtk_reads}} reads. It reads the barcode from a fixed read structure and has no check like the `site` group's; every read here has the remnant, so that check changes no call.
