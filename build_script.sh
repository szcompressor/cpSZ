#!/bin/bash
src_dir=`pwd`

# build FTK for evaluation and speculative compression
mkdir -p ${src_dir}/external
cd ${src_dir}/external
git clone https://github.com/lxAltria/ftk.git
cd ftk
git checkout eval-cpSZ
git reset --hard 78c1e80aa71b52f852eacfdbcea44129e62f3641 
mkdir install
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${src_dir}/external/ftk/install
make -j 4
make install

# build cpSZ
cd ${src_dir}
mkdir build
cd build
cmake ..
make -j 4

