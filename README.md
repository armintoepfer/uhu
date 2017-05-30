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