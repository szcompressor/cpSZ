# cpSZ
Critical point preserving for vector field

# installation
Dependancies: ZSTD, please change the path in CMakeLists.txt
./build_script.sh

# Run compression
./build/bin/sz_compress_cp_preserve_2d_test data/uf.dat data/vf.dat 2400 3600 0.1
decompressed data files are data/uf.dat.out and data/vf.dat.out

# Evaluation
mkdir result && cd result
./../external/ftk/build/bin/ex_cp_2d_bilinear 3600 2400 ../data/uf.dat ../data/vf.dat test
To test 2D linear vector fields, re-build the executables and re-run the compression after changing the corresponding lines (comment line 23 and uncomment line 24) and then evaluate with the following command:
./../external/ftk/build/bin/ex_cp_2d 3600 2400 ../data/uf.dat ../data/vf.dat test 