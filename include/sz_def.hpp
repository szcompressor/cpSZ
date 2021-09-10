#ifndef _sz_def_hpp
#define _sz_def_hpp
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstddef>
#include "sz_huffman.hpp"
#include <vector>

using namespace std;

template <typename T>
union fp_int;

template<>
union fp_int<float>{
	unsigned int ival;
	float fval;
};

template<>
union fp_int<double>{
	unsigned long ival;
	double fval;
};

template <typename T>
inline int 
getExponent(T value);

template<>
inline int 
getExponent<float>(float value){
	fp_int<float> fp;
	fp.fval = value;
	unsigned int ival = fp.ival;
	int exp = (ival & 0x7F800000) >> 23;
	// exp -= 127;
	return exp;	
}

template<>
inline int 
getExponent<double>(double value){
	fp_int<double> fp;
	fp.fval = value;
	unsigned long ival = fp.ival;
	int exp = (int)((ival & 0x7FF0000000000000) >> 52);
	// exp -= 1023;
	return exp;	
}

template <typename T>
inline int 
mantissa_len();

template<>
inline int 
mantissa_len<float>(){
	return 23;
}

template<>
inline int 
mantissa_len<double>(){
	return 52;
}

template <typename T>
inline int 
exp_offset();

template<>
inline int 
exp_offset<float>(){
	return 127;
}

template<>
inline int 
exp_offset<double>(){
	return 1023;
}

struct DSize_2d{
	size_t d1;
	size_t d2;
	size_t num_elements;
	int block_size;
	int max_num_block_elements;
	size_t num_x;
	size_t num_y;
	size_t num_blocks;
	size_t dim0_offset;
	DSize_2d(size_t r1, size_t r2, int bs){
		d1 = r1, d2 = r2;
		num_elements = r1 * r2;
		block_size = bs;
		max_num_block_elements = bs * bs;
		num_x = (r1 - 1) / block_size + 1;
		num_y = (r2 - 1) / block_size + 1;
		num_blocks = num_x * num_y;
		dim0_offset = r2;
	}
};

struct DSize_3d{
	size_t d1;
	size_t d2;
	size_t d3;
	size_t num_elements;
	int block_size;
	int max_num_block_elements;
	size_t num_x;
	size_t num_y;
	size_t num_z;
	size_t num_blocks;
	size_t dim0_offset;
	size_t dim1_offset;
	DSize_3d(size_t r1, size_t r2, size_t r3, int bs){
		d1 = r1, d2 = r2, d3 = r3;
		num_elements = r1 * r2 * r3;
		block_size = bs;
		max_num_block_elements = bs * bs * bs;
		num_x = (r1 - 1) / block_size + 1;
		num_y = (r2 - 1) / block_size + 1;
		num_z = (r3 - 1) / block_size + 1;
		num_blocks = num_x * num_y * num_z;
		dim0_offset = r2 * r3;
		dim1_offset = r3;
	}
};

template<typename T>
struct meanInfo{
	bool use_mean;
	T mean;
	meanInfo(bool use = false, T mean_ = 0){
		use_mean = use;
		mean = mean_;
	}
};

#define MAX(a, b) ((a>b)?(a):(b))
#define MIN(a, b) ((a<b)?(a):(b))
template<class T>
inline T MINF(const T& a, const T& b){
    return (a<b) ? (a) : (b);
}

#define RegCoeffNum2d 3
#define RegCoeffNum3d 4
#define RegErrThreshold 0.1
#define RegCoeffRadius 32768
#define RegCoeffCapacity 65536
#define QuantIntvMeanCapacity 8192
#define QuantIntvSampleDistance 100
#define QuantIntvSampleCapacity 32768
#define QuantIntvAccThreshold 0.999
#define RegThresholdSize2d 8
#define RegThresholdSize3d 4
#define LorenzeNoise2d 0.81
#define LorenzeNoise3d 1.22

#endif

