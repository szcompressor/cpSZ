#include "sz_compress_cp_preserve_2d.hpp"
#include "sz_decompress_cp_preserve_2d.hpp"
#include "sz_lossless.hpp"
#include "utils.hpp"
using namespace std;

int main(int argc, char ** argv){
    size_t num_elements = 0;
    float * U = readfile<float>(argv[1], num_elements);
    float * V = readfile<float>(argv[2], num_elements);
    int r1 = atoi(argv[3]);
    int r2 = atoi(argv[4]);
    double max_eb = atof(argv[5]);

    size_t result_size = 0;
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    cout << "start Compression\n";
    // unsigned char * result = sz_compress_cp_preserve_2d_offline(U, V, r1, r2, result_size, false, max_eb);
    // unsigned char * result = sz_compress_cp_preserve_2d_offline_log(U, V, r1, r2, result_size, false, max_eb);
    // unsigned char * result = sz_compress_cp_preserve_2d_online(U, V, r1, r2, result_size, false, max_eb);
    unsigned char * result = sz_compress_cp_preserve_2d_online_log(U, V, r1, r2, result_size, false, max_eb);
    unsigned char * result_after_lossless = NULL;
    size_t lossless_outsize = sz_lossless_compress(ZSTD_COMPRESSOR, 3, result, result_size, &result_after_lossless);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Compression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    cout << "Compressed size = " << lossless_outsize << ", ratio = " << (2*num_elements*sizeof(float)) * 1.0/lossless_outsize << endl;
    free(result);
    // exit(0);
    err = clock_gettime(CLOCK_REALTIME, &start);
    size_t lossless_output = sz_lossless_decompress(ZSTD_COMPRESSOR, result_after_lossless, lossless_outsize, &result, result_size);
    float * dec_U = NULL;
    float * dec_V = NULL;
    // sz_decompress_cp_preserve_2d_offline<float>(result, r1, r2, dec_U, dec_V);
    // sz_decompress_cp_preserve_2d_offline_log<float>(result, r1, r2, dec_U, dec_V);
    // sz_decompress_cp_preserve_2d_online<float>(result, r1, r2, dec_U, dec_V);
    sz_decompress_cp_preserve_2d_online_log<float>(result, r1, r2, dec_U, dec_V);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Decompression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    free(result_after_lossless);

    verify(U, dec_U, num_elements);
    verify(V, dec_V, num_elements);

    writefile((string(argv[1]) + ".out").c_str(), dec_U, num_elements);
    writefile((string(argv[2]) + ".out").c_str(), dec_V, num_elements);
    free(result);
    free(U);
    free(V);
    free(dec_U);
    free(dec_V);
}