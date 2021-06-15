# cpSZ
Critical point preserving for vector field

# installation
Dependancies: ZSTD, please change the path in CMakeLists.txt
./build_script.sh

# Run compression
./bin/sz_compress_cp_preserve_2d_test uf.dat vf.dat 2400 3600 0.1
decompressed data files are uf.dat.out and vf.dat.out

# Evaluation
./bin/ex_cp_2d 3600 2400 uf.dat vf.dat test
