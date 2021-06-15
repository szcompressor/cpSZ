#ifndef _sz_decompress_pwr_hpp
#define _sz_decompress_pwr_hpp

#include <cstddef>

template<typename T>
T *
sz_decompress_2d_pwr(const unsigned char * compressed, size_t r1, size_t r2);

template<typename T>
T *
sz_decompress_2d_pwr_with_eb(const unsigned char * compressed, const double * log_ebs, size_t r1, size_t r2);

#endif