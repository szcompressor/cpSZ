#ifndef _sz_decompress_cp_preserve_2d_offline_hpp
#define _sz_decompress_cp_preserve_2d_offline_hpp

#include "sz_decompression_utils.hpp"
#include "sz_def.hpp"
#include "sz_prediction.hpp"
#include <vector>

template<typename T>
void
sz_decompress_cp_preserve_2d_offline(const unsigned char * compressed, size_t r1, size_t r2, T *& U, T *& V);

template<typename T>
void
sz_decompress_cp_preserve_2d_offline_log(const unsigned char * compressed, size_t r1, size_t r2, T *& U, T *& V);

template<typename T>
void
sz_decompress_cp_preserve_2d_online(const unsigned char * compressed, size_t r1, size_t r2, T *& U, T *& V);

template<typename T>
void
sz_decompress_cp_preserve_2d_online_fp(const unsigned char * compressed, size_t r1, size_t r2, T *& U, T *& V);

template<typename T>
void
sz_decompress_cp_preserve_2d_online_log(const unsigned char * compressed, size_t r1, size_t r2, T *& U, T *& V);

#endif