#include "sz_compress_cp_preserve_3d.hpp"
#include "sz_decompress_cp_preserve_3d.hpp"
#include "sz_lossless.hpp"
#include "utils.hpp"
#include <cassert>
using namespace std;

int main(int argc, char ** argv){
    size_t num = 0;
    size_t num_tets = 0;
    size_t num_points = 0;
    fprintf(stderr, "reading inputs...\n");
    string file_points = string(argv[1]);
    string file_data = string(argv[2]);
    string file_tets = string(argv[3]);
    double max_eb = atof(argv[4]);
    int * tets = readfile<int>(file_tets.c_str(), num);
    assert(num % 4 == 0);
    num_tets = num / 4;
    printf("num_tets = %d\n", num_tets);
    float * points = readfile<float>(file_points.c_str(), num);
    assert(num % 3 == 0);
    num_points = num / 3;
    printf("num_points = %d\n", num_points);
    float * data = readfile<float>(file_data.c_str(), num);

    size_t result_size = 0;
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    cout << "start Compression\n";
    unsigned char * result = sz_compress_cp_preserve_3d_unstructured(num_points, points, data, num_tets, tets, result_size, max_eb);
    unsigned char * result_after_lossless = NULL;
    size_t lossless_outsize = sz_lossless_compress(ZSTD_COMPRESSOR, 3, result, result_size, &result_after_lossless);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Compression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    cout << "Compressed size = " << lossless_outsize << ", ratio = " << (3*num_points*sizeof(float)) * 1.0/lossless_outsize << endl;
    free(result);
    // exit(0);
    err = clock_gettime(CLOCK_REALTIME, &start);
    size_t lossless_output = sz_lossless_decompress(ZSTD_COMPRESSOR, result_after_lossless, lossless_outsize, &result, result_size);
    float * dec_data = NULL;
    sz_decompress_cp_preserve_3d_unstructured<float>(result, num_points, points, num_tets, tets, dec_data);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Decompression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    free(result_after_lossless);

    verify(data, dec_data, num_points*3);

    writefile((file_data + ".out").c_str(), dec_data, num_points*3);
    free(result);
    free(data);
    free(dec_data);
    free(tets);
    free(points);
}