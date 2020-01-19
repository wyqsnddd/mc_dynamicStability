Dependencies
==

- mc_rtc
- StabiliPlus
- Qhull 7.3.2. It can be installed from source as:
  ```sh
  git clone git@github.com:qhull/qhull.git
  git checkout 2019.1 # git checkout v7.3.2
  cd qhull && mkdir build && cd build
  cmake -DCMAKE_BUILD_TYPE=Release ..
  make -j8
  sudo make install
```
