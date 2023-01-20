#include "sz_cp_preserve_utils.hpp"
#include "sz_compress_3d.hpp"
#include "sz_compress_cp_preserve_2d.hpp"
#include "sz_def.hpp"
#include "sz_compression_utils.hpp"

/*
triangle mesh x0, x1, x2, derive cp-preserving eb for x2 given x0, x1
using SoS method
*/
template<typename T>
double 
derive_cp_abs_eb_sos_online(const T u0, const T u1, const T u2, const T v0, const T v1, const T v2){
	double M0 = u2*v0 - u0*v2;
	double M1 = u1*v2 - u2*v1;
	double M2 = u0*v1 - u1*v0;
	double M = M0 + M1 + M2;
	if(M == 0) return 0;
	bool f1 = (M0 == 0) || (M / M0 >= 1);
	bool f2 = (M1 == 0) || (M / M1 >= 1); 
	bool f3 = (M2 == 0) || (M / M2 >= 1);
	double eb = 0;
	if(f1 && f2 && f3){
		return 0;
	}
	else{
		// test
		// keep sign for the original simplex
		// preserve sign of (M + Ae1 + Be2)
		eb = fabs(M) / (fabs(u1 - u0) + fabs(v0 - v1));
		{
			// keep sign for replacing the first and second vertices
			double cur_eb = MINF(fabs(u1*v2 - u2*v1) / (fabs(u1) + fabs(v1)), fabs(u0*v2 - u2*v0) / (fabs(u0) + fabs(v0)));
			eb = MINF(eb, cur_eb);
		}
	}
	return eb;
}

template<typename T>
unsigned char *
sz_compress_cp_preserve_sos_2d_online(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){
	size_t num_elements = r1 * r2;
	T * decompressed_U = (T *) malloc(num_elements*sizeof(T));
	memcpy(decompressed_U, U, num_elements*sizeof(T));
	T * decompressed_V = (T *) malloc(num_elements*sizeof(T));
	memcpy(decompressed_V, V, num_elements*sizeof(T));
	int * eb_quant_index = (int *) malloc(2*num_elements*sizeof(int));
	int * data_quant_index = (int *) malloc(2*num_elements*sizeof(int));
	int * eb_quant_index_pos = eb_quant_index;
	int * data_quant_index_pos = data_quant_index;
	// next, row by row
	const int base = 2;
	const double log_of_base = log2(base);
	const int capacity = 65536;
	const int intv_radius = (capacity >> 1);
	unpred_vec<T> unpred_data;
	T * U_pos = decompressed_U;
	T * V_pos = decompressed_V;
	// offsets to get six adjacent triangle indices
	// the 7-th rolls back to T0
	/*
			T3	T4
		T2	X 	T5
		T1	T0(T6)
	*/
	const int offsets[7] = {
		-(int)r2, -(int)r2 - 1, -1, (int)r2, (int)r2+1, 1, -(int)r2
	};
	const T x[6][3] = {
		{1, 0, 1},
		{0, 0, 1},
		{0, 1, 1},
		{0, 1, 0},
		{1, 1, 0},
		{1, 0, 0}
	};
	const T y[6][3] = {
		{0, 0, 1},
		{0, 1, 1},
		{0, 1, 0},
		{1, 1, 0},
		{1, 0, 0},
		{1, 0, 1}
	};
	int index_offset[6][2][2];
	for(int i=0; i<6; i++){
		for(int j=0; j<2; j++){
			index_offset[i][j][0] = x[i][j] - x[i][2];
			index_offset[i][j][1] = y[i][j] - y[i][2];
		}
	}
	double threshold = std::numeric_limits<double>::epsilon();
	// conditions_2d cond;
	for(int i=0; i<r1; i++){
		// printf("start %d row\n", i);
		T * cur_U_pos = U_pos;
		T * cur_V_pos = V_pos;
		for(int j=0; j<r2; j++){
			double required_eb = max_pwr_eb;
			// derive eb given six adjacent triangles
			for(int k=0; k<6; k++){
				bool in_mesh = true;
				for(int p=0; p<2; p++){
					// reserved order!
					if(!(in_range(i + index_offset[k][p][1], (int)r1) && in_range(j + index_offset[k][p][0], (int)r2))){
						in_mesh = false;
						break;
					}
				}
				if(in_mesh){
					required_eb = MINF(required_eb, derive_cp_abs_eb_sos_online(cur_U_pos[offsets[k]], cur_U_pos[offsets[k+1]], cur_U_pos[0],
						cur_V_pos[offsets[k]], cur_V_pos[offsets[k+1]], cur_V_pos[0]));
				}
			}
			if(required_eb > 0){
				bool unpred_flag = false;
				T decompressed[2];
				// compress U and V
				for(int k=0; k<2; k++){
					T * cur_data_pos = (k == 0) ? cur_U_pos : cur_V_pos;
					T cur_data = *cur_data_pos;
					double abs_eb = required_eb;
					eb_quant_index_pos[k] = eb_exponential_quantize(abs_eb, base, log_of_base, threshold);
					// eb_quant_index_pos[k] = eb_linear_quantize(abs_eb, 1e-3);
					if(eb_quant_index_pos[k] > 0){
						// get adjacent data and perform Lorenzo
						/*
							d2 X
							d0 d1
						*/
						T d0 = (i && j) ? cur_data_pos[-1 - r2] : 0;
						T d1 = (i) ? cur_data_pos[-r2] : 0;
						T d2 = (j) ? cur_data_pos[-1] : 0;
						T pred = d1 + d2 - d0;
						double diff = cur_data - pred;
						double quant_diff = fabs(diff) / abs_eb + 1;
						if(quant_diff < capacity){
							quant_diff = (diff > 0) ? quant_diff : -quant_diff;
							int quant_index = (int)(quant_diff/2) + intv_radius;
							data_quant_index_pos[k] = quant_index;
							decompressed[k] = pred + 2 * (quant_index - intv_radius) * abs_eb; 
							// check original data
							if(fabs(decompressed[k] - cur_data) >= abs_eb){
								unpred_flag = true;
								break;
							}
						}
						else{
							unpred_flag = true;
							break;
						}
					}
					else unpred_flag = true;
				}
				if(unpred_flag){
					// recover quant index
					*(eb_quant_index_pos ++) = 0;
					*(eb_quant_index_pos ++) = 0;
					*(data_quant_index_pos ++) = intv_radius;
					*(data_quant_index_pos ++) = intv_radius;
					unpred_data.push_back(*cur_U_pos);
					unpred_data.push_back(*cur_V_pos);
				}
				else{
					eb_quant_index_pos += 2;
					data_quant_index_pos += 2;
					// assign decompressed data
					*cur_U_pos = decompressed[0];
					*cur_V_pos = decompressed[1];
				}
			}
			else{
				// record as unpredictable data
				*(eb_quant_index_pos ++) = 0;
				*(eb_quant_index_pos ++) = 0;
				*(data_quant_index_pos ++) = intv_radius;
				*(data_quant_index_pos ++) = intv_radius;
				unpred_data.push_back(*cur_U_pos);
				unpred_data.push_back(*cur_V_pos);
			}
			cur_U_pos ++, cur_V_pos ++;
		}
		U_pos += r2;
		V_pos += r2;
	}
	free(decompressed_U);
	free(decompressed_V);
	printf("offsets eb_q, data_q, unpred: %ld %ld %ld\n", eb_quant_index_pos - eb_quant_index, data_quant_index_pos - data_quant_index, unpred_data.size());
	unsigned char * compressed = (unsigned char *) malloc(2*num_elements*sizeof(T));
	unsigned char * compressed_pos = compressed;
	write_variable_to_dst(compressed_pos, base);
	write_variable_to_dst(compressed_pos, threshold);
	write_variable_to_dst(compressed_pos, intv_radius);
	size_t unpredictable_count = unpred_data.size();
	write_variable_to_dst(compressed_pos, unpredictable_count);
	write_array_to_dst(compressed_pos, (T *)&unpred_data[0], unpredictable_count);	
	Huffman_encode_tree_and_data(2*1024, eb_quant_index, 2*num_elements, compressed_pos);
	free(eb_quant_index);
	Huffman_encode_tree_and_data(2*capacity, data_quant_index, 2*num_elements, compressed_pos);
	printf("pos = %ld\n", compressed_pos - compressed);
	free(data_quant_index);
	compressed_size = compressed_pos - compressed;
	return compressed;	
}

template
unsigned char *
sz_compress_cp_preserve_sos_2d_online(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_sos_2d_online(const double * U, const double * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template<typename T, typename T_fp>
int64_t convert_to_fixed_point(const T * U, const T * V, size_t num_elements, T_fp * U_fp, T_fp * V_fp, T_fp& range, int minbits=8, int maxbits=21){
	double vector_field_resolution = std::numeric_limits<double>::max();
	int64_t vector_field_scaling_factor = 1;
	for (int i=0; i<num_elements; i++){
		double min_val = std::min(fabs(U[i]), fabs(V[i]));
		vector_field_resolution = std::min(vector_field_resolution, min_val);
	}
	int nbits = std::ceil(std::log2(1.0 / vector_field_resolution));
	nbits = std::max(minbits, std::min(nbits, maxbits));
	vector_field_scaling_factor = 1 << nbits;
	std::cerr << "resolution=" << vector_field_resolution 
	<< ", factor=" << vector_field_scaling_factor 
	<< ", nbits=" << nbits << std::endl;
	int64_t max = std::numeric_limits<int64_t>::min();
	int64_t min = std::numeric_limits<int64_t>::max();
	printf("max = %lld, min = %lld\n", max, min);
	for(int i=0; i<num_elements; i++){
		U_fp[i] = U[i] * vector_field_scaling_factor;
		V_fp[i] = V[i] * vector_field_scaling_factor;
		max = std::max(max, U_fp[i]);
		max = std::max(max, V_fp[i]);
		min = std::min(min, U_fp[i]);
		min = std::min(min, V_fp[i]);
	}
	printf("max = %lld, min = %lld\n", max, min);
	range = max - min;
	return vector_field_scaling_factor;
}

// fixed point version
template<typename T_data>
unsigned char *
sz_compress_cp_preserve_sos_2d_online_fp(const T_data * U, const T_data * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){

	using T = int64_t;
	size_t num_elements = r1 * r2;
	T * U_fp = (T *) malloc(num_elements*sizeof(T));
	T * V_fp = (T *) malloc(num_elements*sizeof(T));
	T range = 0;
	T vector_field_scaling_factor = convert_to_fixed_point(U, V, num_elements, U_fp, V_fp, range);
	printf("fixed point range = %lld\n", range);
	int * eb_quant_index = (int *) malloc(2*num_elements*sizeof(int));
	int * data_quant_index = (int *) malloc(2*num_elements*sizeof(int));
	int * eb_quant_index_pos = eb_quant_index;
	int * data_quant_index_pos = data_quant_index;
	// next, row by row
	const int base = 2;
	const double log_of_base = log2(base);
	const int capacity = 65536;
	const int intv_radius = (capacity >> 1);
	T max_eb = range * max_pwr_eb;
	unpred_vec<T> unpred_data;
	T * U_pos = U_fp;
	T * V_pos = V_fp;
	// offsets to get six adjacent triangle indices
	// the 7-th rolls back to T0
	/*
			T3	T4
		T2	X 	T5
		T1	T0(T6)
	*/
	const int offsets[7] = {
		-(int)r2, -(int)r2 - 1, -1, (int)r2, (int)r2+1, 1, -(int)r2
	};
	const int x[6][3] = {
		{1, 0, 1},
		{0, 0, 1},
		{0, 1, 1},
		{0, 1, 0},
		{1, 1, 0},
		{1, 0, 0}
	};
	const int y[6][3] = {
		{0, 0, 1},
		{0, 1, 1},
		{0, 1, 0},
		{1, 1, 0},
		{1, 0, 0},
		{1, 0, 1}
	};
	int index_offset[6][2][2];
	for(int i=0; i<6; i++){
		for(int j=0; j<2; j++){
			index_offset[i][j][0] = x[i][j] - x[i][2];
			index_offset[i][j][1] = y[i][j] - y[i][2];
		}
	}
	T threshold = 1;
	// conditions_2d cond;
	for(int i=0; i<r1; i++){
		// printf("start %d row\n", i);
		T * cur_U_pos = U_pos;
		T * cur_V_pos = V_pos;
		for(int j=0; j<r2; j++){
			T required_eb = max_eb;
			// derive eb given six adjacent triangles
			for(int k=0; k<6; k++){
				bool in_mesh = true;
				for(int p=0; p<2; p++){
					// reserved order!
					if(!(in_range(i + index_offset[k][p][1], (int)r1) && in_range(j + index_offset[k][p][0], (int)r2))){
						in_mesh = false;
						break;
					}
				}
				if(in_mesh){
					required_eb = MINF(required_eb, (T) derive_cp_abs_eb_sos_online(cur_U_pos[offsets[k]], cur_U_pos[offsets[k+1]], cur_U_pos[0],
						cur_V_pos[offsets[k]], cur_V_pos[offsets[k+1]], cur_V_pos[0]));
				}
			}
			if(required_eb > 0){
				bool unpred_flag = false;
				T decompressed[2];
				// compress U and V
				for(int k=0; k<2; k++){
					T * cur_data_pos = (k == 0) ? cur_U_pos : cur_V_pos;
					T cur_data = *cur_data_pos;
					T abs_eb = required_eb;
					eb_quant_index_pos[k] = eb_exponential_quantize(abs_eb, base, log_of_base, threshold);
					// eb_quant_index_pos[k] = eb_linear_quantize(abs_eb, 1e-3);
					if(eb_quant_index_pos[k] > 0){
						// get adjacent data and perform Lorenzo
						/*
							d2 X
							d0 d1
						*/
						T d0 = (i && j) ? cur_data_pos[-1 - r2] : 0;
						T d1 = (i) ? cur_data_pos[-r2] : 0;
						T d2 = (j) ? cur_data_pos[-1] : 0;
						T pred = d1 + d2 - d0;
						T diff = cur_data - pred;
						T quant_diff = std::abs(diff) / abs_eb + 1;
						if(quant_diff < capacity){
							quant_diff = (diff > 0) ? quant_diff : -quant_diff;
							int quant_index = (int)(quant_diff/2) + intv_radius;
							data_quant_index_pos[k] = quant_index;
							decompressed[k] = pred + 2 * (quant_index - intv_radius) * abs_eb; 
							// check original data
							if(std::abs(decompressed[k] - cur_data) >= abs_eb){
								unpred_flag = true;
								break;
							}
						}
						else{
							unpred_flag = true;
							break;
						}
					}
					else unpred_flag = true;
				}
				if(unpred_flag){
					// recover quant index
					*(eb_quant_index_pos ++) = 0;
					*(eb_quant_index_pos ++) = 0;
					*(data_quant_index_pos ++) = intv_radius;
					*(data_quant_index_pos ++) = intv_radius;
					unpred_data.push_back(*cur_U_pos);
					unpred_data.push_back(*cur_V_pos);
				}
				else{
					eb_quant_index_pos += 2;
					data_quant_index_pos += 2;
					// assign decompressed data
					*cur_U_pos = decompressed[0];
					*cur_V_pos = decompressed[1];
				}
			}
			else{
				// record as unpredictable data
				*(eb_quant_index_pos ++) = 0;
				*(eb_quant_index_pos ++) = 0;
				*(data_quant_index_pos ++) = intv_radius;
				*(data_quant_index_pos ++) = intv_radius;
				unpred_data.push_back(*cur_U_pos);
				unpred_data.push_back(*cur_V_pos);
			}
			cur_U_pos ++, cur_V_pos ++;
		}
		U_pos += r2;
		V_pos += r2;
	}
	free(U_fp);
	free(V_fp);
	printf("offsets eb_q, data_q, unpred: %ld %ld %ld\n", eb_quant_index_pos - eb_quant_index, data_quant_index_pos - data_quant_index, unpred_data.size());
	unsigned char * compressed = (unsigned char *) malloc(2*num_elements*sizeof(T));
	unsigned char * compressed_pos = compressed;
	write_variable_to_dst(compressed_pos, vector_field_scaling_factor);
	write_variable_to_dst(compressed_pos, base);
	write_variable_to_dst(compressed_pos, threshold);
	write_variable_to_dst(compressed_pos, intv_radius);
	size_t unpredictable_count = unpred_data.size();
	write_variable_to_dst(compressed_pos, unpredictable_count);
	write_array_to_dst(compressed_pos, (T *)&unpred_data[0], unpredictable_count);	
	Huffman_encode_tree_and_data(2*1024, eb_quant_index, 2*num_elements, compressed_pos);
	free(eb_quant_index);
	Huffman_encode_tree_and_data(2*capacity, data_quant_index, 2*num_elements, compressed_pos);
	printf("pos = %ld\n", compressed_pos - compressed);
	free(data_quant_index);
	compressed_size = compressed_pos - compressed;
	return compressed;	

}

template
unsigned char *
sz_compress_cp_preserve_sos_2d_online_fp(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_sos_2d_online_fp(const double * U, const double * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);
