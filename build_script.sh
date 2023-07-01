#!/bin/bash
src_dir=`pwd`

# build FTK for evaluation and speculative compression
mkdir -p ${src_dir}/external
cd ${src_dir}/external
git clone https://github.com/hguo/ftk.git
cd ftk
#git reset --hard 5af454b9471c0f617d3d170446b42051c0d0a4d6
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

