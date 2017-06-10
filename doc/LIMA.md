<h1 align="center">
    lima - CCS Barcode Demultiplexer
</h1>

<p align="center">
  <img src="img/lima.png" alt="Logo of Lima" width="150px"/>
</p>

## TOC
* [Scope](#scope)
* [Barcode score](#barcode-score)
* [Defaults](#defaults)
* [Barcoding modes](#barcoding-modes)
* [Alignment options](#alignment-options)

## Scope
Demultiplexes CCS reads with insane speed, vectorized alignment and parallelized
processing. In- and output are BAM. Barcode sequences get clipped and `bq` and `bc` tags
added, just like bam2bam. Barcodes do not necessarily have to be in the correct
direction.

## Output
*Lima* generates three output files, all starting with the BAM input file name
prefix.

### BAM
The first file `prefix.demux.bam` contains clipped subreads, annotated with
barcode tags, that passed filters.

### Report
Second file is `prefix.demux.report`, a tsv file about each read, unfiltered.
Example:

    $ head prefix.demux.report | column -t
    ZMW                                BcLeft  BcRight  ScoreLeft  ScoreRight  Score  ClipLeft  ClipRight
    m54011_170105_093642/30867881/ccs  0       50       84         59          71     14        2223
    m54011_170105_093642/30867884/ccs  36      14       78         100         89     15        2222
    m54011_170105_093642/30867886/ccs  3       36       47         100         73     15        2214
    m54011_170105_093642/30867887/ccs  10      32       100        100         100    15        2217

### Summary
Third file is `prefix.demux.summary`, showing how many reads have been filtered:

    Above length and score threshold : 979
    Below length and score threshold : 2
    Below length threshold           : 5
    Below score threshold            : 21

## Barcode score
The barcode score is normalized and between 0 and 100, whereas 0 is no hit and
100 perfect match. The provided score is the mean of both normalizated
barcode scores:

    score = (left_barcode_ssw_score + right_barcode_ssw_score) / (2 * barcode_length * match_score)

## Defaults
 - Reads with length below 50 bp after demultiplexing are omitted.
   Adjusted with `--min-length`
 - Reads with barcode score below 50 are omitted.
   Adjust with `--min-score`
 - For each barcode, we align it to a subsequence of the begin and end of
   the CCS read. The length of the subsequence is `barcode_length * multiplier`,
   which can be adjusted with `--window-size-mult`.

## Barcoding modes
Currently we support following modes:

### Symmetric
For symmetric, please use

    --mode symmetric

If your barcodes are not in correct direction, please try

    --mode symmetric --try-rc

### Tailed
For tailed, please use

    --mode symmetric --try-rc

### Asymmetric
For asymmetric, we try every barcode as given and reverse complement for the
left and right side of the ccs read separately, find the best matching barcode,
and report them together. Option `--try-rc` is implictly activated, please use

    --mode asymmetric

## Alignment options

    -A,--match-score       Score for a sequence match. [2]
    -B,--mismatch-penalty  Penalty for a mismatch. [2]
    -O,--gap-open-penalty  Gap open penalties for deletions and insertions. [3]
    -e,--gap-ext-penalty   Gap extension penalties for deletions and insertions. [1]