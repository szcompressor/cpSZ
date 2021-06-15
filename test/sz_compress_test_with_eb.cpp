#include "sz_compress_pwr.hpp"
#include "sz_decompress_pwr.hpp"
#include "sz_compress_3d.hpp"
#include "sz_decompress_3d.hpp"
#include "sz_lossless.hpp"
#include "sz_compression_utils.hpp"
#include "utils.hpp"
#include <limits>
using namespace std;

int main(int argc, char ** argv){
    size_t num_elements = 0;
    double * data = readfile<double>(argv[1], num_elements);
    int r1 = atoi(argv[2]);
    int r2 = atoi(argv[3]);
    int r3 = atoi(argv[4]);
    int log_base = atoi(argv[5]);
    // double max = data[0];
    // double min = data[0];
    // for(int i=1; i<num_elements; i++){
    //     if(max < data[i]) max = data[i];
    //     if(min > data[i]) min = data[i];
    // }
    // cout << "value range = " << max - min << endl;
    // cout << "precision = " << eb*(max - min) << endl;
    size_t result_size = 0;
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    double * precisions = (double *) malloc(num_elements*sizeof(double));
    // for(int i=0; i<num_elements; i++){
    //     precisions[i] = (double)rand()/RAND_MAX;
    //     precisions[i] = MAX(1e-5, precisions[i]);
    // }
    // {
    size_t tmp_num = 0;
    double * tmp = readfile<double>("/Users/xin/github/ftk_xin/ftk/build/ebs.dat", tmp_num);
    for(int i=0; i<num_elements; i++){
        precisions[i] = tmp[i] * fabs(data[i]);
    }
    // quantize precisions
    size_t inds_lossless_outsize = 0;
    {
        int inc = log_base;
        double threshold = std::numeric_limits<float>::epsilon();
        int * inds = (int *) malloc(num_elements*sizeof(int));
        for(int i=0; i<num_elements; i++){
            if(precisions[i] < threshold){
                precisions[i] = 0;
                inds[i] = 0;
            }
            else{
                int id = log2(precisions[i] / threshold)/log2(inc);
                if(id > 255) id = 255;
                inds[i] = id;
                precisions[i] = pow(inc, inds[i]) * threshold;
            }
        }
        unsigned char * tmp_inds = (unsigned char *) malloc(num_elements*sizeof(int));
        unsigned char * tmp_inds_pos = tmp_inds;
        Huffman_encode_tree_and_data(2*256, inds, num_elements, tmp_inds_pos);
        free(inds);
        unsigned char * tmp_inds_lossless = NULL;
        inds_lossless_outsize = sz_lossless_compress(ZSTD_COMPRESSOR, 3, tmp_inds, tmp_inds_pos - tmp_inds, &tmp_inds_lossless);
        cout << "Error bounds compressed size = " << inds_lossless_outsize << endl;
        free(tmp_inds_lossless);
        free(tmp_inds);
    }
    // free(tmp);
    // }
    // double * log_precisions = NULL;
    unsigned char * result = sz_compress_2d_with_eb(data, precisions, r1*r2, r3, result_size);
    // unsigned char * result = sz_compress_2d_pwr_with_eb(data, precisions, &log_precisions, r1*r2, r3, result_size);
    unsigned char * result_after_lossless = NULL;
    size_t lossless_outsize = sz_lossless_compress(ZSTD_COMPRESSOR, 3, result, result_size, &result_after_lossless);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Compression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    cout << "Compressed size = " << lossless_outsize << ", ratio = " << (num_elements*sizeof(double)) * 1.0/lossless_outsize << endl;
    cout << "Final compression ratio = " << (num_elements*sizeof(double)) * 1.0/(lossless_outsize + inds_lossless_outsize) << endl;
    free(result);
    err = clock_gettime(CLOCK_REALTIME, &start);
    size_t lossless_output = sz_lossless_decompress(ZSTD_COMPRESSOR, result_after_lossless, lossless_outsize, &result, result_size);
    double * dec_data = sz_decompress_2d_with_eb<double>(result, precisions, r1*r2, r3);
    // double * dec_data = sz_decompress_2d_pwr_with_eb<double>(result, log_precisions, r1*r2, r3);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Decompression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    free(result_after_lossless);
    // writefile("dec_data.dat", dec_data, num_elements);
    for(int i=0; i<num_elements; i++){
        if(precisions[i] > 0){
            if(fabs((dec_data[i] - data[i]) / data[i]) > tmp[i]){
                // cout << i << " " << precisions[i] << " " << log_precisions[i] << ": " << data[i] << " " << dec_data[i] << " " << fabs((dec_data[i] - data[i]) / data[i]) << endl;
                cout << i << " " << precisions[i] << " " << tmp[i] << " " << tmp[i]*data[i] << ": " << data[i] << " " << dec_data[i] << " " << fabs((dec_data[i] - data[i]) / data[i]) << endl;
                exit(0);
            } 
        }
    }
    free(tmp);
    string fn = string(argv[1]) + ".out";
    // {
    //     size_t size = num_elements;
    //     float * uf = (float *) malloc(size*sizeof(float));
    //     for(int i=0; i<size; i++){
    //         uf[i] = dec_data[i];
    //     }
    //     writefile(fn.c_str(), uf, size);
    //     free(uf);
    // }
    writefile(fn.c_str(), dec_data, num_elements);
    verify(data, dec_data, num_elements);
    free(precisions);
    // free(log_precisions);
    free(result);
    free(data);
    free(dec_data);
}