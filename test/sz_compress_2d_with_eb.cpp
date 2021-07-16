#include "sz_compress_3d.hpp"
#include "sz_decompress_3d.hpp"
#include "sz_lossless.hpp"
#include "utils.hpp"
#include "sz_cp_preserve_utils.hpp"
#include "sz_compression_utils.hpp"
#include <cassert>
#include <vector>
using namespace std;

size_t quantize_and_compress_eb(double * eb, size_t num_elements, double base_eb=1e-7){
    int * quant_ind = (int *) malloc(num_elements * sizeof(int));
    const int base = 2;
    const double log_of_base = log2(base);
    for(int i=0; i<num_elements; i++){
        quant_ind[i] = eb_exponential_quantize(eb[i], base, log_of_base, base_eb);    
    }
    unsigned char * tmp = (unsigned char *) malloc(num_elements * sizeof(int));
    unsigned char * tmp_pos = tmp;
    Huffman_encode_tree_and_data(2*1024, quant_ind, num_elements, tmp_pos);
    size_t compressed_eb_size = tmp_pos - tmp;
    unsigned char * tmp2 = (unsigned char *) malloc(num_elements * sizeof(int));
    size_t lossless_outsize = sz_lossless_compress(ZSTD_COMPRESSOR, 3, tmp, compressed_eb_size, &tmp2);
    free(tmp);
    free(quant_ind);
    free(tmp2);
    return lossless_outsize;
}

int main(int argc, char ** argv){
    size_t num_elements = 0;
    float * data = readfile<float>(argv[1], num_elements);
    size_t num_elements_eb = 0;
    double * eb = readfile<double>(argv[2], num_elements_eb);
    assert(num_elements == num_elements_eb);
    int r1 = atoi(argv[3]);
    int r2 = atoi(argv[4]);
    bool quantize_eb = atoi(argv[5]);
    double base_eb = 1e-9;
    if(argc > 6) base_eb = atof(argv[6]);

    size_t eb_size = 0;
    if(quantize_eb){
        eb_size = quantize_and_compress_eb(eb, num_elements, base_eb);
        std::cout << "Compressed eb size = " << eb_size << endl;
    }

    size_t result_size = 0;
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    unsigned char * result =  sz_compress_2d_with_eb(data, eb, r1, r2, result_size);
    unsigned char * result_after_lossless = NULL;
    size_t lossless_outsize = sz_lossless_compress(ZSTD_COMPRESSOR, 3, result, result_size, &result_after_lossless);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Compression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    cout << "Data compressed size = " << lossless_outsize << ", ratio = " << (num_elements*sizeof(float)) * 1.0/lossless_outsize << endl;
    if(quantize_eb){
        cout << "Overall compressed size = " << lossless_outsize + eb_size << ", ratio = " << (num_elements*sizeof(float)) * 1.0/(lossless_outsize + eb_size) << endl;
    }
    free(result);
    err = clock_gettime(CLOCK_REALTIME, &start);
    size_t lossless_output = sz_lossless_decompress(ZSTD_COMPRESSOR, result_after_lossless, lossless_outsize, &result, result_size);
    float * dec_data = sz_decompress_2d_with_eb<float>(result, eb, r1, r2);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Decompression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    free(result_after_lossless);
    // writefile("dec_data.dat", dec_data, num_elements);
    for(int i=0; i<num_elements; i++){
        if(fabs(data[i] - dec_data[i]) > eb[i]){
            cerr << "Error bound is not respected in " << i << "-th element: " << data[i] << " " << dec_data[i] << std::endl;
            exit(0);
        }
    }
    verify(data, dec_data, num_elements);
    free(result);
    free(data);
    free(dec_data);
}