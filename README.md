<p align="center">
  <img src="doc/img/uhu.png" alt="uhu logos" width="150px"/>
</p>
<h1 align="center">UHU</h1>
<p align="center">Sandbox Tools for PacBio Reads</p>

***

## INSTALL

  ```sh
  git clone https://github.com/armintoepfer/uhu.git && cd uhu
  git submodule update --init
  git submodule foreach git pull origin develop
  mkdir build && cd build
  cmake -GNinja -DCMAKE_INSTALL_PREFIX=~/bin .. && ninja
  ninja install
  ```

## TOOLS

### zmw_stats
Extracts the following metrics for each ZMW to standard out:
 - ZMW ID
 - Number of counts for each CX tag between 0 - 15
 - Forward barcode
 - Barcode Score
 - Number of Subreads
 - Subread Mean Length
 - Subread Median Length
 - Subread Length Standard Deviation

### demux_ccs
Demultiplexes CCS reads. Output is BAM. Only symmetric mode is supported.
Reads below 50 bp after demultiplexing are omitted. Barcode sequences get
removed and `bq` and `bc` tags added.