#include "sz_decompress_3d.hpp"
#include "sz_decompress_cp_preserve_2d.hpp"
#include "sz_decompress_block_processing.hpp"
#include <limits>
#include <unordered_set>

template<typename T>
void
sz_decompress_cp_preserve_3d_online_log(const unsigned char * compressed, size_t r1, size_t r2, size_t r3, T *& U, T *& V, T *& W){
	if(U) free(U);
	if(V) free(V);
	if(W) free(W);
	size_t num_elements = r1 * r2 * r3;
	const unsigned char * compressed_pos = compressed;
	int base = 0;
	read_variable_from_src(compressed_pos, base);
	printf("base = %d\n", base);
	int intv_radius = 0;
	read_variable_from_src(compressed_pos, intv_radius);
	size_t sign_map_size = (num_elements - 1)/8 + 1;
	unsigned char * sign_map_u = convertByteArray2IntArray_fast_1b_sz(num_elements, compressed_pos, sign_map_size);	
	unsigned char * sign_map_v = convertByteArray2IntArray_fast_1b_sz(num_elements, compressed_pos, sign_map_size);	
	unsigned char * sign_map_w = convertByteArray2IntArray_fast_1b_sz(num_elements, compressed_pos, sign_map_size);	
	const int capacity = (intv_radius << 1);
	size_t unpred_data_count = 0;
	read_variable_from_src(compressed_pos, unpred_data_count);
	const T * eb_zero_data = (T *) compressed_pos;
	const T * eb_zero_data_pos = eb_zero_data;
	compressed_pos += unpred_data_count*sizeof(T);
	size_t eb_quant_num = 0;
	read_variable_from_src(compressed_pos, eb_quant_num);
	int * eb_quant_index = Huffman_decode_tree_and_data(2*256, eb_quant_num, compressed_pos);
	size_t data_quant_num = 0;
	read_variable_from_src(compressed_pos, data_quant_num);
	int * data_quant_index = Huffman_decode_tree_and_data(2*capacity, data_quant_num, compressed_pos);
	U = (T *) malloc(num_elements*sizeof(T));
	V = (T *) malloc(num_elements*sizeof(T));
	W = (T *) malloc(num_elements*sizeof(T));
	T * U_pos = U;
	T * V_pos = V;
	T * W_pos = W;
	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	int * eb_quant_index_pos = eb_quant_index;
	int * data_quant_index_pos = data_quant_index;
	const double threshold=std::numeric_limits<float>::epsilon();
	double log_of_base = log2(base);
	int eb_quant_index_max = (int) (log2(1.0 / threshold)/log_of_base) + 1;
	std::unordered_set<int> unpred_data_indices;
	for(int i=0; i<r1; i++){
		for(int j=0; j<r2; j++){
			for(int k=0; k<r3; k++){
				// printf("%ld %ld %ld\n", i, j, k);
				T * data_pos[3] = {U_pos, V_pos, W_pos};
				int index = i*dim0_offset + j*dim1_offset + k;
				// get eb
				if(*eb_quant_index_pos == 0 || *eb_quant_index_pos == eb_quant_index_max){
					unpred_data_indices.insert(index);
					for(int p=0; p<3; p++){
						T cur_data = *(eb_zero_data_pos ++);
						*(data_pos[p]) = (cur_data == 0) ? -100 : log2f(fabs(cur_data));
					}
					eb_quant_index_pos ++;
				}
				else{
					double eb = (*eb_quant_index_pos == 0) ? 0 : pow(base, *eb_quant_index_pos) * threshold;
					eb_quant_index_pos ++;
					for(int p=0; p<3; p++){
						T * cur_log_data_pos = data_pos[p];					
						T d0 = (i && j && k) ? cur_log_data_pos[- dim0_offset - dim1_offset - 1] : 0;
						T d1 = (i && j) ? cur_log_data_pos[- dim0_offset - dim1_offset] : 0;
						T d2 = (i && k) ? cur_log_data_pos[- dim0_offset - 1] : 0;
						T d3 = (i) ? cur_log_data_pos[- dim0_offset] : 0;
						T d4 = (j && k) ? cur_log_data_pos[- dim1_offset - 1] : 0;
						T d5 = (j) ? cur_log_data_pos[- dim1_offset] : 0;
						T d6 = (k) ? cur_log_data_pos[- 1] : 0;
						T pred = d0 + d3 + d5 + d6 - d1 - d2 - d4;
						*cur_log_data_pos = pred + 2 * (data_quant_index_pos[p] - intv_radius) * eb;
					}
					data_quant_index_pos += 3;
				}
				U_pos ++;
				V_pos ++;
				W_pos ++;
			}
		}
	}
	printf("recover data done\n");
	eb_zero_data_pos = eb_zero_data;
	for(int i=0; i<num_elements; i++){
		if(unpred_data_indices.count(i)){
			U[i] = *(eb_zero_data_pos++);
			V[i] = *(eb_zero_data_pos++);
			W[i] = *(eb_zero_data_pos++);
		}
		else{
			if(U[i] < -99) U[i] = 0;
			else U[i] = sign_map_u[i] ? exp2(U[i]) : -exp2(U[i]);
			if(V[i] < -99) V[i] = 0;
			else V[i] = sign_map_v[i] ? exp2(V[i]) : -exp2(V[i]);
			if(W[i] < -99) W[i] = 0;
			else W[i] = sign_map_w[i] ? exp2(W[i]) : -exp2(W[i]);
		}
	}
	free(sign_map_u);
	free(sign_map_v);
	free(sign_map_w);
	free(eb_quant_index);
	free(data_quant_index);
}

template
void
sz_decompress_cp_preserve_3d_online_log<float>(const unsigned char * compressed, size_t r1, size_t r2, size_t r3, float *& U, float *& V, float *& W);

template
void
sz_decompress_cp_preserve_3d_online_log<double>(const unsigned char * compressed, size_t r1, size_t r2, size_t r3, double *& U, double *& V, double *& W);

template<typename T>
void
sz_decompress_cp_preserve_3d_unstructured(const unsigned char * compressed, int n, const T * points, int m, const int * tets_ind, T *& data){

	if(data) free(data);
	const unsigned char * compressed_pos = compressed;
	int base = 0;
	read_variable_from_src(compressed_pos, base);
	printf("base = %d\n", base);
	int intv_radius = 0;
	read_variable_from_src(compressed_pos, intv_radius);
	size_t sign_map_size = (3*n - 1)/8 + 1;
	unsigned char * sign_map = convertByteArray2IntArray_fast_1b_sz(3*n, compressed_pos, sign_map_size);	
	const int capacity = (intv_radius << 1);
	size_t unpred_data_count = 0;
	read_variable_from_src(compressed_pos, unpred_data_count);
	const T * eb_zero_data = (T *) compressed_pos;
	const T * eb_zero_data_pos = eb_zero_data;
	compressed_pos += unpred_data_count*sizeof(T);
	size_t eb_quant_num = 0;
	read_variable_from_src(compressed_pos, eb_quant_num);
	int * eb_quant_index = Huffman_decode_tree_and_data(2*256, eb_quant_num, compressed_pos);
	size_t data_quant_num = 0;
	read_variable_from_src(compressed_pos, data_quant_num);
	int * data_quant_index = Huffman_decode_tree_and_data(2*capacity, data_quant_num, compressed_pos);
	data = (T *) malloc(3*n*sizeof(T));
	T * data_pos = data;
	int * eb_quant_index_pos = eb_quant_index;
	int * data_quant_index_pos = data_quant_index;
	const double threshold=std::numeric_limits<float>::epsilon();
	double log_of_base = log2(base);
	int eb_quant_index_max = (int) (log2(1.0 / threshold)/log_of_base) + 1;
	std::unordered_set<int> unpred_data_indices;
	for(int i=0; i<n; i++){
		// printf("%ld %ld %ld\n", i, j, k);
		// get eb
		if(*eb_quant_index_pos == 0 || *eb_quant_index_pos == eb_quant_index_max){
			unpred_data_indices.insert(i);
			for(int p=0; p<3; p++){
				T cur_data = *(eb_zero_data_pos ++);
				data_pos[p] = (cur_data == 0) ? -100 : log2f(fabs(cur_data));
			}
			eb_quant_index_pos ++;
		}
		else{
			double eb = (*eb_quant_index_pos == 0) ? 0 : pow(base, *eb_quant_index_pos) * threshold;
			eb_quant_index_pos ++;
			for(int p=0; p<3; p++){
				T * cur_log_data_pos = data_pos + p;					
				T pred = (i) ? cur_log_data_pos[- 3] : 0;
				*cur_log_data_pos = pred + 2 * (data_quant_index_pos[p] - intv_radius) * eb;
			}
			data_quant_index_pos += 3;
		}
		data_pos += 3;
	}
	printf("recover data done\n");
	eb_zero_data_pos = eb_zero_data;
	unsigned char * sign_pos = sign_map;
	for(int i=0; i<n; i++){
		if(unpred_data_indices.count(i)){
			for(int p=0; p<3; p++){
				data[3*i + p] = *(eb_zero_data_pos++);
			}
			sign_pos += 3;
		}
		else{
			for(int p=0; p<3; p++){
				if(data[3*i + p] < -99) data[3*i + p] = 0;
				else data[3*i + p] = *(sign_pos ++) ? exp2(data[3*i + p]) : -exp2(data[3*i + p]);
			}
		}
	}
	free(sign_map);
	free(eb_quant_index);
	free(data_quant_index);
}
template
void
sz_decompress_cp_preserve_3d_unstructured<float>(const unsigned char * compressed, int n, const float * points, int m, const int * tets_ind, float *& data);

template
void
sz_decompress_cp_preserve_3d_unstructured<double>(const unsigned char * compressed, int n, const double * points, int m, const int * tets_ind, double *& data);
