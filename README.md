<p align="center">
  <img src="doc/img/uhu.png" alt="uhu logos" width="150px"/>
</p>
<h1 align="center">UHU</h1>
<p align="center">Sandbox for PacBio Tools</p>

***

## TOOLS

- [Demultiplex PacBio Reads `lima`](doc/LIMA.md)
- [Extract ZMW Stats `zmw_stats`](doc/ZMW_STATS.md)

## INSTALL

Building from scratch requires system-wide installed boost (>=1.58.0),
cmake (3.2), and a c++14 compiler (>=gcc-6, clang). `ninja` or
`make` is used to build (ninja is faster than make).

  ```sh
  git clone https://github.com/armintoepfer/uhu.git && cd uhu
  git submodule update --init
  git submodule foreach git pull origin develop
  mkdir build && cd build
  cmake -GNinja -DCMAKE_INSTALL_PREFIX=~/ .. && ninja
  ninja install
  ```

Always start with a fresh build directory after pulling new commits.