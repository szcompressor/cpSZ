#!/bin/bash
src_dir=`pwd`

# build cpSZ
cd ${src_dir}
mkdir build
cd build
cmake ..
make -j 4

# build FTK for evaluation
mkdir -p ${src_dir}/external
cd ${src_dir}/external
git clone https://github.com/lxAltria/ftk.git
cd ftk
git checkout eval-cpSZ
git reset --hard 78c1e80aa71b52f852eacfdbcea44129e62f3641 
mkdir build
cd build
cmake ..
make -j 4
