#include "sz_compress_cp_preserve_3d.hpp"
#include "sz_decompress_cp_preserve_3d.hpp"
#include "sz_lossless.hpp"
#include "utils.hpp"
using namespace std;

int main(int argc, char ** argv){
    size_t num_elements = 0;
    float * U = readfile<float>(argv[1], num_elements);
    float * V = readfile<float>(argv[2], num_elements);
    float * W = readfile<float>(argv[3], num_elements);
    int r1 = atoi(argv[4]);
    int r2 = atoi(argv[5]);
    int r3 = atoi(argv[6]);
    double max_eb = atof(argv[7]);
    // cout << U[r2 + 3] << " " << U[3*r2 + 1] << endl;
    // transpose_2d(U, r1, r2);
    // cout << U[r2 + 3] << " " << U[3*r2 + 1] << endl;
    // transpose_2d(V, r1, r2);

    size_t result_size = 0;
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    cout << "start Compression\n";
    // unsigned char * result =  sz_compress_cp_preserve_3d_offline_log(U, V, W, r1, r2, r3, result_size, false, max_eb);
    unsigned char * result =  sz_compress_cp_preserve_3d_online_log(U, V, W, r1, r2, r3, result_size, false, max_eb);
    // exit(0);
    unsigned char * result_after_lossless = NULL;
    size_t lossless_outsize = sz_lossless_compress(ZSTD_COMPRESSOR, 3, result, result_size, &result_after_lossless);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Compression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    cout << "Compressed size = " << lossless_outsize << ", ratio = " << (3*num_elements*sizeof(float)) * 1.0/lossless_outsize << endl;
    free(result);
    // exit(0);
    err = clock_gettime(CLOCK_REALTIME, &start);
    size_t lossless_output = sz_lossless_decompress(ZSTD_COMPRESSOR, result_after_lossless, lossless_outsize, &result, result_size);
    float * dec_U = NULL;
    float * dec_V = NULL;
    float * dec_W = NULL;
    sz_decompress_cp_preserve_3d_online_log<float>(result, r1, r2, r3, dec_U, dec_V, dec_W);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Decompression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    free(result_after_lossless);

    verify(U, dec_U, num_elements);
    verify(V, dec_V, num_elements);
    verify(W, dec_W, num_elements);

    // transpose_2d(dec_U, r1, r2);
    // transpose_2d(dec_V, r1, r2);
    writefile((string(argv[1]) + ".out").c_str(), dec_U, num_elements);
    writefile((string(argv[2]) + ".out").c_str(), dec_V, num_elements);
    writefile((string(argv[3]) + ".out").c_str(), dec_W, num_elements);
    free(result);
    free(U);
    free(V);
    free(W);
    free(dec_U);
    free(dec_V);
    free(dec_W);
}