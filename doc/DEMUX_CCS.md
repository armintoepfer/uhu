<h1 align="center">
    demux_ccs
</h1>

## Scope
Demultiplexes CCS reads with insane speed, vectorized alignment and parallelized
processing. Output is BAM. Barcode sequences get clipped and `bq` and `bc` tags
added, just like bam2bam. Barcodes do not necessarily have to be in the correct
direction.

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
 - For each barcode, we align it to a subsequence of the beginning and end of
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
and report them together; please use

    --mode asymmetric