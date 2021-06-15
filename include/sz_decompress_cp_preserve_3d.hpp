#ifndef _sz_decompress_cp_preserve_3d_offline_hpp
#define _sz_decompress_cp_preserve_3d_offline_hpp

#include "sz_decompression_utils.hpp"
#include "sz_def.hpp"
#include "sz_prediction.hpp"
#include <vector>

template<typename T>
void
sz_decompress_cp_preserve_3d_online_log(const unsigned char * compressed, size_t r1, size_t r2, size_t r3, T *& U, T *& V, T *& W);

template<typename T>
void
sz_decompress_cp_preserve_3d_unstructured(const unsigned char * compressed, int n, const T * points, int m, const int * tets_ind, T *& data);

#endif