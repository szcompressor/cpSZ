#!/bin/bash
src_dir=`pwd`
# build cpSZ
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
git commit 11830bb5b4469ed1e7c94c3c3235a2ebfbd46145 
mkdir build
cd build
cmake ..
make -j 4
