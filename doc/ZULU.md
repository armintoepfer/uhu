<h1 align="center">
    zulu - Barcode validator
</h1>

<p align="center">
  <img src="img/zulu.png" alt="Logo of Zulu" width="150px"/>
</p>

## Scope
If you have a known mapping of Barcode IDs to unique references, `zulu` is
your weapon of choice. `Zulu` computes various metrics to quality control
your barcodes:
 - Positive Predictive Value
 - Number of ZMWs mapped
 - Optimal barcode score given a PPV

## Workflow
Align your barcoded reads against the known templates. Each template name must
be unique. Example:

    $ cat refs.fasta
    >reference.A
    ACAGGACGTACG...
    >reference.B
    ATATTATTTTAT...
    >reference.C
    CAGCAGACTTGC...

    $ blasr movie.barcoded.bam refs.fasta --bam --out aligned.bam

Prepare a list of barcode ID to template. Barcodes start at index 0.
Your example list of assignments:

    0=reference.A,1=reference.B,2=reference.C

Call *zulu* with your mapping and a minimal aligned read length that you trust:

    zulu aligned.bam --min-length XXX --mapping "0=reference.A,1=reference.B,2=reference.C"

## Output
Following an exemplary output with explanations:

    #Subreads input        : 9545  <- All reads in the input file
    #Subreads BC & >1500bp : 7670  <- Reads that are barcoded and above min-length

    #ZMWs input            : 980   <- All ZMWs in the input file
    #ZMWs BC & >1500bp     : 906   <- ZMWs with barcoded and above min-length subreads

    PPV                    : 0.91  <- Positive predictive value

## Advanced Options
### Bam2bam Input
If you demultiplexed using bam2bam and your library was tailed or asymmetric use `--tailed`

### Compute PPV per ZMW
Take one and only one subread per ZMW `--zmw`

### Compute False Negative rate
Provide the number of expectet barcode IDs, e.g. 384, `--num-barcodes 384`

### Compute Minimal Barcode Score to Reach PPV
Provide the PPV of interest, e.g. 0.99, `--min-ppv 0.99`