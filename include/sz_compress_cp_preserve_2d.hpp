#ifndef _sz_compress_cp_preserve_2d_hpp
#define _sz_compress_cp_preserve_2d_hpp

#include <cstddef>

#define DEFAULT_EB 0.1

template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_offline(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose=false, double max_pwr_eb=0.1);

template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_offline_log(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose=false, double max_pwr_eb=0.1);

template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_online(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose=false, double max_pwr_eb=0.1);

// template<typename T>
// unsigned char *
// sz_compress_cp_preserve_sos_2d_online(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose=false, double max_pwr_eb=0.1);

template<typename T>
unsigned char *
sz_compress_cp_preserve_sos_2d_online_fp(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose=false, double max_pwr_eb=0.1);

template<typename T>
unsigned char *
sz_compress_cp_preserve_sos_2d_online_fp_spec_eb(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose=false, double max_pwr_eb=0.1);

template<typename T>
unsigned char *
sz_compress_cp_preserve_sos_2d_online_fp_spec_exec_fn(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose=false, double max_pwr_eb=0.1, double max_factor=1);

template<typename T>
unsigned char *
sz_compress_cp_preserve_sos_2d_online_fp_spec_exec_all(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose=false, double max_pwr_eb=0.1, double max_factor=1);

template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_online_log(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose=false, double max_pwr_eb=0.1);

template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_bilinear_online_log(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose=false, double max_pwr_eb=0.1);

#endif