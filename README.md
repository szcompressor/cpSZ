# cpSZ
Critical point preserving for vector field

# installation
sh build_script.sh <br>
Dependancies: ZSTD

# Run compression
./build/bin/sz_compress_cp_preserve_2d_test data/uf.dat data/vf.dat 2400 3600 0.1 <br>
Decompressed data files are data/uf.dat.out and data/vf.dat.out

# Evaluation
./external/ftk/build/bin/ex_cp_2d_bilinear 3600 2400 ./data/uf.dat ./data/vf.dat test <br>
To test 2D linear vector fields, re-build the executables and re-run the compression after changing the corresponding lines (comment line 23 and uncomment line 24) and then evaluate with the following command: <br>
./external/ftk/build/bin/ex_cp_2d 3600 2400 ./data/uf.dat ./data/vf.dat test 