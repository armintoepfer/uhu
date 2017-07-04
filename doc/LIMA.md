<h1 align="center">
    lima - PacBio Barcode Demultiplexer
</h1>

<p align="center">
  <img src="img/lima.png" alt="Logo of Lima" width="150px"/>
</p>

## TOC
* [Scope](#scope)
* [Workflow](#workflow)
* [Barcode score](#barcode-score)
* [Defaults](#defaults)
* [Missing Features](#missing-features)
* [FAQ](#faq)
  * [How fast is fast?](#how-fast-is-fast)
  * [Is there a way to show the progress?](#is-there-a-way-to-show-the-progress)
  * [Can I set the number of threads?](#can-i-set-the-number-of-threads)
  * [How can I easily plot the score distributions?](#how-can-i-easily-plot-the-score-distributions)
  * [Can I split my data by barcode?](#can-i-split-my-data-by-barcode)
  * [Why are asymmetric hits reported in --symmetric mode?](#why-are-asymmetric-hits-reported-in---symmetric-mode)
  * [Why are symmetric hits reported in the default asymmetric mode?](#why-are-symmetric-hits-reported-in-the-default-asymmetric-mode)
  * [How do barcode indices correspond to the input sequences?](#how-do-barcode-indices-correspond-to-the-input-sequences)
  * [I used the tailed library prep, what options to choose?](#i-used-the-tailed-library-prep,-what-options-to-choose)
  * [What is different to bam2bam?](#what-is-different-to-bam2bam)

## Scope
*Lima* offers following features:
 * Demultiplex PacBio reads with insane speed, vectorized alignment and parallelized processing
 * Both, raw subreads and ccs reads can be processed
 * In- and output are BAM
 * Barcode sequences get clipped and `bq` and `bc` tags added, just like bam2bam
 * Barcodes do not necessarily have to be in the correct direction
 * Output can be split by barcode
 * No scraps.bam needed

## Execution
Run on raw subread data:

    lima movie.subreads.bam barcodes.fasta

Run on CCS data:

    lima --css movie.ccs.bam barcodes.fasta

## Workflow

<img src="img/barcode.png" width="1000px">

*Lima* processes input reads grouped by ZMW.
Each target barcode region, left and right, is processed individually.
For a particular target barcode region, every barcode sequence gets
aligned as given and as reverse-complement and per barcode scores are summed;
the best scoring barcode is chosen using the score sums.

This procedure corresponds to the *asymmetric* library prep.
If only identical barcode pairs are of interest, *symmetric*, please use
`--symmetric`.

## Output
Both *lima* tools generate four output files, all starting with the BAM input
file name prefix.

### BAM
The first file `prefix.demux.bam` contains clipped records, annotated with
barcode tags, that passed filters and respects `--symmetric`.

### Report
Second file is `prefix.demux.report`, a tsv file about each read, unfiltered.
An individual score with `-1` indicates that a leading or trailing adapter is
missing. This is irrelevant for CCS reads.

    $ head prefix.demux.report | column -t
    ZMW      IndexLeft  IndexRight  MeanScoreLeft  MeanScoreRight  MeanScore  ClipsLeft    ClipsRight           ScoresLeft    ScoresRight
    4391559  2          2           73             100             87         0,14,15,14  1558,2097,2183,2113  -1,56,82,82   100,100,100,-1
    4457329  2          2           65             85              75         0,15,18     2772,2174,2402       -1,54,76      87,82,-1
    4522785  3          3           86             87              87         0,15,15,14  2016,2176,2198,2119  -1,100,76,82  73,100,89,-1

### Summary
Third file is `prefix.demux.summary`, shows how many ZMWs have been filtered,
how ZMWs many are *symmetric*/*asymmetric*, and how many reads have been filtered.

    ZMWs above length and score threshold : 1127
    ZMWs below length and score threshold : 0
    ZMWs below length threshold           : 0
    ZMWs below score threshold            : 2025

    ZMWs symmetric                        : 1013
    ZMWs asymmetric                       : 114

    Reads above length                    : 7596
    Reads below length                    : 9

### Counts
Fourth file is `prefix.demux.counts`, a tsv file, shows the counts for each
observed barcode pair; only those barcode that passed filters are counted.
Example:

    $ cat prefix.demux.counts | column -t
    IndexLeft  IndexRight  Counts
    0          0           100
    0          1           56
    1          0           5112
    1          2           846

## Barcode score
The barcode score is normalized between 0 and 100, whereas 0 is no hit and
100 perfect match. The provided mean score is the mean of both normalizated
barcode scores.

## Defaults
 - Reads with length below 50 bp after demultiplexing are omitted.
   Adjusted with `--min-length`.
 - ZMWs with barcode score below 50 are omitted.
   Adjust with `--min-score`
 - For each barcode, we align it to a subsequence of the begin and end of
   the CCS read. The length of the subsequence is `barcode_length * multiplier`,
   which can be adjusted with `--window-size-mult`.
 - Alignment options
    -A,--match-score       Score for a sequence match.
    -B,--mismatch-penalty  Penalty for a mismatch.
    -O,--gap-open-penalty  Gap open penalties for deletions and insertions.
    -e,--gap-ext-penalty   Gap extension penalties for deletions and insertions.

## Missing Features
 * bam2bam-like BAM barcode header line
 * verification
 * validation

## FAQ
### How fast is fast?
Example: 200 barcodes, asymmetric mode (try each barcode forward and
reverse-complement), 300,000 CCS reads. On my 2014 iMac with 4 cores + HT:

    503.57s user 11.74s system 725% cpu 1:11.01 total

Those 8:09 minutes translate into 0.233 milliseconds per ZMW,
1.16 microseconds per barcode for both sides aligning forward and reverse-complement,
and 291 nanoseconds per alignment. This includes IO.

### Is there a way to show the progress?
No. Please run `wc -l prefix.demux.report` to get the number of processed ZMWs.

### How can I easily plot the score distributions?
Use `R`. Example:

    r = read.table("~/Downloads/ccs.demux.report", header = TRUE)
    par(mfrow=c(3,1))
    hist(r$ScoreLeft,breaks=0:100,xlab="",ylab="Counts",main="Left Barcode Score")
    hist(r$ScoreRight,breaks=0:100,xlab="",ylab="Counts",main="Right Barcode Score")
    hist(r$MeanScore,breaks=0:100,xlab="Barcode Score",ylab="Counts",main="Combined Average Barcode Score")

<img src="img/score_hist.png" width="1000px">

### Can I split my data by barcode?
You can either iterate of the `prefix.demux.bam` file N times or use
`--split-bam`. Each barcode has its own BAM file called
`prefix.leftIdx-rightIdx.demux.bam`, e.g., `prefix.0-0.demux.bam`.
This mode consumes more memory, as output cannot be streamed.

### Why are asymmetric hits reported in --symmetric mode?
*Lima* tries all barcode combinations and `--symmetric` only filters BAM output.
Sequences flanked by *asymmetric* barcodes are still reported, but are not
written to BAM. By not enforcing only *symmetric* barcode pairs, *lima* gains
higher PPV, as your sample might be contaminated and contains unwanted
barcode pairs; instead of enforcing one *symmetric* pair, *lima* rather
filters such sequences. Every *symmetric* library contains few *asymmetric*
templates. If many *asymmetric* templates are called, your library preparation
might be bad.

### Why are symmetric hits reported in the default asymmetric mode?
Even if your sample is labeled *asymmetric*, *symmetric* hits are simply
sequences flanked by the same barcode.

### How do barcode indices correspond to the input sequences?
Input barcode sequences are tagged with an incrementing counter. The first
sequence is barcode `0` and the last barcode `numBarcodes - 1`.

### I used the tailed library prep, what options to choose?
Use `--symmetric`.

### What is different to bam2bam?
 * CCS read support
 * Barcodes of every adapter gets scored for raw subreads
 * Do not enforce symmetric barcode pairing, which increases PPV
 * For asymmetric barcodes, `lima` reports the identified order, instead of
   descending sorting
 * Call barcodes per read and not per adapter
 * Open-source and can be compiled on your local Mac or Linux machine
 * Faster
 * Nice reports for QC