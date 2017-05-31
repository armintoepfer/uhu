<p align="center">
  <img src="doc/img/uhu.png" alt="uhu logos" width="150px"/>
</p>
<h1 align="center">UHU</h1>
<p align="center">Sandbox Tools for PacBio Reads</p>

***

## INSTALL

  ```sh
  git clone ssh://git@bitbucket.nanofluidics.com:7999/~atoepfer/uhu.git && cd uhu
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
Demultiplexes CCS reads. Only symmetric mode supported.