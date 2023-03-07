#ifndef _sz_compression_utils_hpp
#define _sz_compression_utils_hpp

#include "sz_def.hpp"

template <typename T>
inline T
out_of_range_data_encode(const T data, int residue_len, unsigned char *& sign_pos, int *& type_pos, unsigned char *& dst, int& pos){
	*(sign_pos ++) = (data>0);
	*(type_pos ++) = getExponent(data);
	fp_int<T> fp;
	fp.fval = data;
	int discard_len = mantissa_len<T>() - residue_len;
	fp.ival = (fp.ival >> discard_len) << discard_len;
	T decompressed = fp.fval;
	fp.ival &= 0x000FFFFFFFFFFFFF;
	fp.ival >>= discard_len;
	while(residue_len){
		int byte_rest_len = 8 - pos;
		if(residue_len >= byte_rest_len){
			*dst = (*dst) | (0xFF & (fp.ival >> (residue_len - byte_rest_len)));
			residue_len -= byte_rest_len;
			dst ++;
			pos = 0;
		}
		else{
			*dst = (*dst) | (((0xFF >> (8 - residue_len)) & fp.ival) << (byte_rest_len - residue_len));
			pos += residue_len;
			break;
		}
	}
	return decompressed;
}

template <typename T>
inline void
write_variable_to_dst(unsigned char *& dst, const T& var){
    memcpy(dst, &var, sizeof(T));
    dst += sizeof(T);
}

template <typename T>
inline void
write_array_to_dst(unsigned char *& dst, const T * array, size_t length){
    memcpy(dst, array, length*sizeof(T));
    dst += length*sizeof(T);
}

// quantize with decompression data (Lorenzo)
template<typename T>
inline int
quantize(float pred, T cur_data, double precision, int capacity, int intv_radius, T *& unpredictable_data_pos, T * decompressed){
	double diff = cur_data - pred;
	double quant_diff = fabs(diff) / precision + 1;
	if(quant_diff < capacity){
		quant_diff = (diff > 0) ? quant_diff : -quant_diff;
		int quant_index = (int)(quant_diff/2) + intv_radius;
		T decompressed_data = pred + 2 * (quant_index - intv_radius) * precision; 
		*decompressed = decompressed_data;
		if(fabs(decompressed_data - cur_data) <= precision) return quant_index;
 	}
 	*decompressed = cur_data;
 	*(unpredictable_data_pos++) = cur_data;
 	return 0;
}

// return quantization index, no decompression data (regression)
template<typename T>
inline int
quantize(float pred, T cur_data, double precision, int capacity, int intv_radius, T *& unpredictable_data_pos){
	double diff = cur_data - pred;
	double quant_diff = fabs(diff) / precision + 1;
	if(quant_diff < capacity){
		quant_diff = (diff > 0) ? quant_diff : -quant_diff;
		int quant_index = (int)(quant_diff/2) + intv_radius;
		T decompressed_data = pred + 2 * (quant_index - intv_radius) * precision; 
		if(fabs(decompressed_data - cur_data) <= precision) return quant_index;
 	}
 	*(unpredictable_data_pos++) = cur_data;
 	return 0;
}

inline void
compress_regression_coefficient_2d(const double * precisions, float * reg_params_pos, int * reg_params_type_pos, float *& reg_unpredictable_data_pos){
	float * prev_reg_params = reg_params_pos - RegCoeffNum2d;
	for(int i=0; i<RegCoeffNum2d; i++){
		*(reg_params_type_pos ++) = quantize(*prev_reg_params, *reg_params_pos, precisions[i], RegCoeffCapacity, RegCoeffRadius, reg_unpredictable_data_pos, reg_params_pos);
		prev_reg_params ++, reg_params_pos ++; 
	}
}

inline void
compress_regression_coefficient_3d(const double * precisions, float * reg_params_pos, int * reg_params_type_pos, float *& reg_unpredictable_data_pos){
	float * prev_reg_params = reg_params_pos - RegCoeffNum3d;
	for(int i=0; i<RegCoeffNum3d; i++){
		*(reg_params_type_pos ++) = quantize(*prev_reg_params, *reg_params_pos, precisions[i], RegCoeffCapacity, RegCoeffRadius, reg_unpredictable_data_pos, reg_params_pos);
		prev_reg_params ++, reg_params_pos ++; 
	}
}

void
encode_regression_coefficients_2d(const int * reg_params_type, const float * reg_unpredictable_data, size_t reg_count, size_t reg_unpredictable_count, unsigned char *& compressed_pos);

void
encode_regression_coefficients(const int * reg_params_type, const float * reg_unpredictable_data, size_t reg_count, size_t reg_unpredictable_count, unsigned char *& compressed_pos);

// copied from conf.c
unsigned int 
round_up_power_of_2(unsigned int base);

// modified from TypeManager.c
// change return value and increment byteArray
void 
convertIntArray2ByteArray_fast_1b_to_result_sz(const unsigned char* intArray, size_t intArrayLength, unsigned char *& compressed_pos);

// for test meta-compressor. sz implementation can remove this header.
HuffmanTree *
build_Huffman_tree(const int * type, size_t num_elements, size_t state_num);

void
Huffman_encode_tree_and_data(size_t state_num, const int * type, size_t num_elements, unsigned char*& compressed_pos);

// variation with speculative compression on derived eb
template<typename T>
T relax_eb(T eb, T factor){
	return eb * factor;
}

template<typename T>
T relax_eb(T eb){
	return eb * 2;
}

template<typename T>
T restrict_eb(T eb){
	return eb / 2;
}


#endif
