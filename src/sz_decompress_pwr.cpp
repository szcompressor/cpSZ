#include "sz_decompress_pwr.hpp"
#include "sz_decompress_3d.hpp"
#include "sz_decompression_utils.hpp"
#include <cmath>
#include <cstring>

template<typename T>
T *
sz_decompress_2d_pwr(const unsigned char * compressed, size_t r1, size_t r2){
	size_t num_elements = r1 * r2;
	size_t sign_size = (num_elements-1)/8 + 1;
	const unsigned char * compressed_pos = compressed;
	T threshold = *((T *) compressed_pos);
	compressed_pos += sizeof(T);
	unsigned char * signs = convertByteArray2IntArray_fast_1b_sz(num_elements, compressed_pos, sign_size);	
	T * dec_data = sz_decompress_2d<T>(compressed_pos, r1, r2);
	for(size_t i=0; i<num_elements; i++){
		if(dec_data[i] < threshold) dec_data[i] = 0;
		else dec_data[i] = exp2(dec_data[i]);
		if(signs[i]) dec_data[i] = -dec_data[i];
	}
	free(signs);
	return dec_data;
}
template
float * 
sz_decompress_2d_pwr<float>(const unsigned char * compressed, size_t r1, size_t r2);

template<typename T>
T *
sz_decompress_2d_pwr_with_eb(const unsigned char * compressed, const double * log_ebs, size_t r1, size_t r2){
	size_t num_elements = r1 * r2;
	size_t sign_size = (num_elements-1)/8 + 1;
	const unsigned char * compressed_pos = compressed;
	T threshold = *((T *) compressed_pos);
	compressed_pos += sizeof(T);
	unsigned char * signs = convertByteArray2IntArray_fast_1b_sz(num_elements, compressed_pos, sign_size);	
	printf("threshold = %.4g\nBegin decompress log data\n", threshold);
	T * dec_data = sz_decompress_2d_with_eb<T>(compressed_pos, log_ebs, r1, r2);
	for(size_t i=0; i<num_elements; i++){
		if(dec_data[i] == threshold) dec_data[i] = 0;
		else dec_data[i] = exp2((double)dec_data[i]);
		if(signs[i]) dec_data[i] = -dec_data[i];
	}
	free(signs);
	return dec_data;
}
template
float * 
sz_decompress_2d_pwr_with_eb<float>(const unsigned char * compressed, const double * log_ebs, size_t r1, size_t r2);
