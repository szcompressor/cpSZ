#include "sz_cp_preserve_utils.hpp"
#include "sz_compress_cp_preserve_3d.hpp"
#include "sz_def.hpp"
#include "sz_compression_utils.hpp"
#include <unordered_map>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/numeric/fixed_point.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/print.hh>

// offsets to get 24 adjacent simplex indices
// x -> z, fastest -> slowest dimensions
// current data would always be the last index, i.e. x[i][3]
static const int coordinates[24][4][3] = {
	// offset = 0, 0, 0
	{
		{0, 0, 1},
		{0, 1, 1},
		{1, 1, 1},
		{0, 0, 0}
	},
	{
		{0, 1, 0},
		{0, 1, 1},
		{1, 1, 1},
		{0, 0, 0}
	},
	{
		{0, 0, 1},
		{1, 0, 1},
		{1, 1, 1},
		{0, 0, 0}
	},
	{
		{1, 0, 0},
		{1, 0, 1},
		{1, 1, 1},
		{0, 0, 0}
	},
	{
		{0, 1, 0},
		{1, 1, 0},
		{1, 1, 1},
		{0, 0, 0}
	},
	{
		{1, 0, 0},
		{1, 1, 0},
		{1, 1, 1},
		{0, 0, 0}
	},
	// offset = -1, 0, 0
	{
		{0, 0, 0},
		{1, 0, 1},
		{1, 1, 1},
		{1, 0, 0}
	},
	{
		{0, 0, 0},
		{1, 1, 0},
		{1, 1, 1},
		{1, 0, 0}
	},
	// offset = 0, -1, 0
	{
		{0, 0, 0},
		{0, 1, 1},
		{1, 1, 1},
		{0, 1, 0}
	},
	{
		{0, 0, 0},
		{1, 1, 0},
		{1, 1, 1},
		{0, 1, 0}
	},
	// offset = -1, -1, 0
	{
		{0, 0, 0},
		{0, 1, 0},
		{1, 1, 1},
		{1, 1, 0}
	},
	{
		{0, 0, 0},
		{1, 0, 0},
		{1, 1, 1},
		{1, 1, 0}
	},
	// offset = 0, 0, -1
	{
		{0, 0, 0},
		{0, 1, 1},
		{1, 1, 1},
		{0, 0, 1}
	},
	{
		{0, 0, 0},
		{1, 0, 1},
		{1, 1, 1},
		{0, 0, 1}
	},
	// offset = -1, 0, -1
	{
		{0, 0, 0},
		{0, 0, 1},
		{1, 1, 1},
		{1, 0, 1}
	},
	{
		{0, 0, 0},
		{1, 0, 0},
		{1, 1, 1},
		{1, 0, 1}
	},
	// offset = 0, -1, -1
	{
		{0, 0, 0},
		{0, 0, 1},
		{1, 1, 1},
		{0, 1, 1}
	},
	{
		{0, 0, 0},
		{0, 1, 0},
		{1, 1, 1},
		{0, 1, 1}
	},
	// offset = -1, -1, -1
	{
		{0, 0, 0},
		{0, 0, 1},
		{0, 1, 1},
		{1, 1, 1}
	},
	{
		{0, 0, 0},
		{0, 1, 0},
		{0, 1, 1},
		{1, 1, 1}
	},
	{
		{0, 0, 0},
		{0, 0, 1},
		{1, 0, 1},
		{1, 1, 1}
	},
	{
		{0, 0, 0},
		{1, 0, 0},
		{1, 0, 1},
		{1, 1, 1}
	},
	{
		{0, 0, 0},
		{0, 1, 0},
		{1, 1, 0},
		{1, 1, 1}
	},
	{
		{0, 0, 0},
		{1, 0, 0},
		{1, 1, 0},
		{1, 1, 1}
	}
};

template<typename T_fp>
static int 
check_cp(T_fp vf[4][3], int indices[4]){
	// robust critical point test
	bool succ = ftk::robust_critical_point_in_simplex3(vf, indices);
	if (!succ) return -1;
	return 1;
}

template<typename T_fp_acc, typename T_fp>
static inline void 
update_index(T_fp_acc vf[4][3], int indices[4], int local_id, int global_id, const T_fp * U, const T_fp * V, const T_fp * W){
	indices[local_id] = global_id;
	vf[local_id][0] = U[global_id];
	vf[local_id][1] = V[global_id];
	vf[local_id][2] = W[global_id];
}

template<typename T_fp>
static vector<bool> 
compute_cp(const T_fp * U_fp, const T_fp * V_fp, const T_fp * W_fp, int r1, int r2, int r3){
	// check cp for all cells
	vector<bool> cp_exist(6*(r1-1)*(r2-1)*(r3-1), 0);
	// vector<bool> cp_exist(6*r1*r2*r3, 0);
	ptrdiff_t dim0_offset = r2*r3;
	ptrdiff_t dim1_offset = r3;
	ptrdiff_t cell_dim0_offset = (r2-1)*(r3-1);
	ptrdiff_t cell_dim1_offset = r3-1;
	int indices[4] = {0};
	__int128 vf[4][3] = {0};
	for(int i=0; i<r1-1; i++){
		for(int j=0; j<r2-1; j++){
			for(int k=0; k<r3-1; k++){
				bool verbose = false;
				// order (reserved, z->x):
				ptrdiff_t cell_offset = 6*(i*cell_dim0_offset + j*cell_dim1_offset + k);
				// ptrdiff_t tmp_offset = 6*(i*dim0_offset + j*dim1_offset + k);
				// if(tmp_offset == 4353132) verbose = true;
				// if(verbose){
				// 	std::cout << i << " " << j << " " << k << std::endl;
				// }
				// (ftk-0) 000, 001, 011, 111
				update_index(vf, indices, 0, i*dim0_offset + j*dim1_offset + k, U_fp, V_fp, W_fp);
				update_index(vf, indices, 1, (i+1)*dim0_offset + j*dim1_offset + k, U_fp, V_fp, W_fp);
				update_index(vf, indices, 2, (i+1)*dim0_offset + (j+1)*dim1_offset + k, U_fp, V_fp, W_fp);
				update_index(vf, indices, 3, (i+1)*dim0_offset + (j+1)*dim1_offset + (k+1), U_fp, V_fp, W_fp);
				cp_exist[cell_offset] = (check_cp(vf, indices) == 1);
				// if(verbose){
				// 	auto offset = cell_offset;
				// 	std::cout << "cell id = " << offset << ", cp = " << +cp_exist[offset] << std::endl;
				// 	std::cout << "indices: ";
				// 	for(int i=0; i<4; i++){
				// 		std::cout << indices[i] << " "; 
				// 	}
				// 	std::cout << std::endl;
				// 	T_fp tmp[4][3];
				// 	for(int i=0; i<4; i++){
				// 		for(int j=0; j<3; j++){
				// 			tmp[i][j] = vf[i][j];
				// 		}
				// 	}
				// 	ftk::print4x3("M:", tmp);
				// }				
				// (ftk-2) 000, 010, 011, 111
				update_index(vf, indices, 1, i*dim0_offset + (j+1)*dim1_offset + k, U_fp, V_fp, W_fp);
				cp_exist[cell_offset + 1] = (check_cp(vf, indices) == 1);
				// (ftk-1) 000, 001, 101, 111
				update_index(vf, indices, 1, (i+1)*dim0_offset + j*dim1_offset + k, U_fp, V_fp, W_fp);
				update_index(vf, indices, 2, (i+1)*dim0_offset + j*dim1_offset + k+1, U_fp, V_fp, W_fp);
				cp_exist[cell_offset + 2] = (check_cp(vf, indices) == 1);
				// (ftk-4) 000, 100, 101, 111
				update_index(vf, indices, 1, i*dim0_offset + j*dim1_offset + k+1, U_fp, V_fp, W_fp);
				cp_exist[cell_offset + 3] = (check_cp(vf, indices) == 1);
				// (ftk-3) 000, 010, 110, 111
				update_index(vf, indices, 1, i*dim0_offset + (j+1)*dim1_offset + k, U_fp, V_fp, W_fp);
				update_index(vf, indices, 2, i*dim0_offset + (j+1)*dim1_offset + k+1, U_fp, V_fp, W_fp);
				cp_exist[cell_offset + 4] = (check_cp(vf, indices) == 1);
				// if(verbose){
				// 	auto offset = cell_offset + 4;
				// 	std::cout << "cell id = " << offset << ", cp =" << +cp_exist[offset] << std::endl;
				// 	std::cout << "indices: ";
				// 	for(int i=0; i<4; i++){
				// 		std::cout << indices[i] << " "; 
				// 	}
				// 	std::cout << std::endl;
				// 	T_fp tmp[4][3];
				// 	for(int i=0; i<4; i++){
				// 		for(int j=0; j<3; j++){
				// 			tmp[i][j] = vf[i][j];
				// 		}
				// 	}
				// 	ftk::print4x3("M:", tmp);
				// }				
				// (ftk-5) 000, 100, 110, 111
				update_index(vf, indices, 1, i*dim0_offset + j*dim1_offset + k+1, U_fp, V_fp, W_fp);
				cp_exist[cell_offset + 5] = (check_cp(vf, indices) == 1);
			}
		}
	}
	return cp_exist;	
}

template <typename T> 
static bool 
same_direction(T u0, T u1, T u2, T u3) {
    int sgn0 = sgn(u0);
    if(sgn0 == 0) return false;
    if((sgn0 == sgn(u1)) && (sgn0 == sgn(u2)) && (sgn0 == sgn(u3))) return true;
    return false;
}

template<typename T>
static inline T
det_2_by_2(const T u0, const T u1, const T v0, const T v1){
	return u0*v1 - u1*v0;
}

template<typename T>
static inline T
det_3_by_3(const T u0, const T u1, const T u2, const T v0, const T v1, const T v2, const T w0, const T w1, const T w2){
	return u0*v1*w2 + u1*v2*w0 + u2*v0*w1 - u0*v2*w1 - u1*v0*w2 -u2*v1*w0;
}

template<typename T>
static inline T
abs(const T u){
	return (u>0) ? u : -u;
}
/*
Tet x0, x1, x2, x3, derive cp-preserving eb for x3 given x0, x1, x2
using SoS method
*/
template<typename T_acc, typename T>
static T 
derive_cp_abs_eb_sos_online(const T u0, const T u1, const T u2, const T u3, const T v0, const T v1, const T v2, const T v3, const T w0, const T w1, const T w2, const T w3){
	T_acc M0 = - det_3_by_3(u1, u2, u3, v1, v2, v3, w1, w2, w3);
	T_acc M1 = + det_3_by_3(u0, u2, u3, v0, v2, v3, w0, w2, w3);
	T_acc M2 = - det_3_by_3(u0, u1, u3, v0, v1, v3, w0, w1, w3);
	T_acc M3 = + det_3_by_3(u0, u1, u2, v0, v1, v2, w0, w1, w2);
	T_acc M = M0 + M1 + M2 + M3;
	if(M == 0) return 0;
	T same_eb = 0;
	if(same_direction(u0, u1, u2, u3)){			
		same_eb = MAX(same_eb, std::abs(u3));
	}
	if(same_direction(v0, v1, v2, v3)){			
		same_eb = MAX(same_eb, std::abs(v3));
	}
	if(same_direction(w0, w1, w2, w3)){			
		same_eb = MAX(same_eb, std::abs(w3));
	}
	if(same_eb != 0) return same_eb;
	// keep sign for the original simplex
	T one = 1;
	T denominator = abs(det_3_by_3<T_acc>(v0, v1, v2, w0, w1, w2, one, one, one)) + abs(det_3_by_3<T_acc>(u0, u1, u2, w0, w1, w2, one, one, one)) + abs(det_3_by_3<T_acc>(u0, u1, u2, v0, v1, v2, one, one, one)); 
	T eb = abs(M) / denominator;
	{
		// keep sign for replacing the three other vertices
		T cur_eb_0 = abs(M0)/(std::abs(det_2_by_2(v1, v2, w1, w2)) + std::abs(det_2_by_2(u1, u2, w1, w2)) + std::abs(det_2_by_2(v1, v2, w1, w2)));
		T cur_eb_1 = abs(M1)/(std::abs(det_2_by_2(v0, v2, w0, w2)) + std::abs(det_2_by_2(u0, u2, w0, w2)) + std::abs(det_2_by_2(v0, v2, w0, w2)));
		T cur_eb_2 = abs(M2)/(std::abs(det_2_by_2(v0, v1, w0, w1)) + std::abs(det_2_by_2(u0, u1, w0, w1)) + std::abs(det_2_by_2(v0, v1, w0, w1)));
		eb = MINF(MINF(cur_eb_0, cur_eb_1), MINF(cur_eb_2, eb));
	}
	return eb;
}

template<typename T, typename T_fp>
static int64_t 
convert_to_fixed_point(const T * U, const T * V, const T * W, size_t num_elements, T_fp * U_fp, T_fp * V_fp, T_fp * W_fp, T_fp& range, int minbits=8, int maxbits=23){
	double vector_field_resolution = std::numeric_limits<double>::max();
	int64_t vector_field_scaling_factor = 1;
	for (int i=0; i<num_elements; i++){
		double min_val = std::min(std::min(fabs(U[i]), fabs(V[i])), fabs(W[i]));
		vector_field_resolution = std::min(vector_field_resolution, min_val);
	}
	int nbits = maxbits;
	if(vector_field_resolution) nbits = std::ceil(std::log2(1.0 / vector_field_resolution));
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
		W_fp[i] = W[i] * vector_field_scaling_factor;
		max = std::max(max, U_fp[i]);
		max = std::max(max, V_fp[i]);
		max = std::max(max, W_fp[i]);
		min = std::min(min, U_fp[i]);
		min = std::min(min, V_fp[i]);
		min = std::min(min, W_fp[i]);
	}
	printf("max = %lld, min = %lld\n", max, min);
	range = max - min;
	// {
	// 	// test 
	// 	// (id==20588) || (id==21100) || (id==21101) || (id==283245)
	// 	int tmp[4] = {20588, 21100, 21101, 283245};
	// 	for(int i=0; i<4; i++){
	// 		auto id = tmp[i];
	// 		std::cout << U_fp[id] << " " << V_fp[id] << " " << W_fp[id] << std::endl;
	// 		T U = U_fp[id] * (T)1.0 / vector_field_scaling_factor;
	// 		T V = V_fp[id] * (T)1.0 / vector_field_scaling_factor;
	// 		T W = W_fp[id] * (T)1.0 / vector_field_scaling_factor;
	// 		T_fp U_ = U * vector_field_scaling_factor;
	// 		T_fp V_ = V * vector_field_scaling_factor;
	// 		T_fp W_ = W * vector_field_scaling_factor;
	// 		std::cout << U_ << " " << V_ << " " << W_ << std::endl;
	// 	}
	// }
	return vector_field_scaling_factor;
}

template<typename T_data>
unsigned char *
sz_compress_cp_preserve_sos_3d_online_fp(const T_data * U, const T_data * V, const T_data * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb){
	std::cout << "sz_compress_cp_preserve_sos_3d_online_fp" << std::endl;
	using T = int64_t;
	size_t num_elements = r1 * r2 * r3;
	T * U_fp = (T *) malloc(num_elements*sizeof(T));
	T * V_fp = (T *) malloc(num_elements*sizeof(T));
	T * W_fp = (T *) malloc(num_elements*sizeof(T));
	T range = 0;
	T vector_field_scaling_factor = convert_to_fixed_point(U, V, W, num_elements, U_fp, V_fp, W_fp, range);
	printf("fixed point range = %lld\n", range);
	int * eb_quant_index = (int *) malloc(num_elements*sizeof(int));
	int * data_quant_index = (int *) malloc(3*num_elements*sizeof(int));
	int * eb_quant_index_pos = eb_quant_index;
	int * data_quant_index_pos = data_quant_index;
	// next, row by row
	const int base = 2;
	const double log_of_base = log2(base);
	const int capacity = 65536;
	const int intv_radius = (capacity >> 1);
	T max_eb = range * max_pwr_eb;
	unpred_vec<T_data> unpred_data;
	ptrdiff_t dim0_offset = r2 * r3;
	ptrdiff_t dim1_offset = r3;
	ptrdiff_t cell_dim0_offset = (r2-1) * (r3-1);
	ptrdiff_t cell_dim1_offset = r3-1;
	ptrdiff_t simplex_offset[24];
	{
		ptrdiff_t * simplex_offset_pos = simplex_offset;
		ptrdiff_t base = 0;
		// offset = 0, 0, 0
		for(int i=0; i<6; i++){
			*(simplex_offset_pos++) = i;
		}
		// offset = -1, 0, 0
		base = -6;
		*(simplex_offset_pos++) = base + 3;
		*(simplex_offset_pos++) = base + 5;
		// offset = 0, -1, 0
		base = -6*cell_dim1_offset;
		*(simplex_offset_pos++) = base + 1;
		*(simplex_offset_pos++) = base + 4;
		// offset = -1, -1, 0
		base = -6 - 6*cell_dim1_offset;
		*(simplex_offset_pos++) = base + 4;
		*(simplex_offset_pos++) = base + 5;
		// offset = 0, 0, -1
		base = -6*cell_dim0_offset;
		*(simplex_offset_pos++) = base + 0;
		*(simplex_offset_pos++) = base + 2;
		// offset = -1, 0, -1
		base = -6*cell_dim0_offset - 6;
		*(simplex_offset_pos++) = base + 2;
		*(simplex_offset_pos++) = base + 3;
		// offset = 0, -1, -1
		base = -6*cell_dim1_offset - 6*cell_dim0_offset;
		*(simplex_offset_pos++) = base + 0;
		*(simplex_offset_pos++) = base + 1;
		// offset = -1, -1, -1
		base = -6*cell_dim0_offset - 6*cell_dim1_offset - 6;
		for(int i=0; i<6; i++){
			*(simplex_offset_pos++) = base + i;
		}
	}
	int index_offset[24][3][3];
	for(int i=0; i<24; i++){
		for(int j=0; j<3; j++){
			for(int k=0; k<3; k++){
				index_offset[i][j][k] = coordinates[i][j][k] - coordinates[i][3][k];
			}
		}
	}
	ptrdiff_t offset[24][3];
	for(int i=0; i<24; i++){
		for(int x=0; x<3; x++){
			// offset[i][x] = (coordinates[i][x][0] - coordinates[i][3][0]) * dim0_offset + (coordinates[i][x][1] - coordinates[i][3][1]) * dim1_offset + (coordinates[i][x][2] - coordinates[i][3][2]);
			offset[i][x] = (coordinates[i][x][0] - coordinates[i][3][0]) + (coordinates[i][x][1] - coordinates[i][3][1]) * dim1_offset + (coordinates[i][x][2] - coordinates[i][3][2]) * dim0_offset;
		}
	}
	T * cur_U_pos = U_fp;
	T * cur_V_pos = V_fp;
	T * cur_W_pos = W_fp;
	T threshold = 1;
	// check cp for all cells
	std::cout << "start cp checking\n";
	vector<bool> cp_exist = compute_cp(U_fp, V_fp, W_fp, r1, r2, r3);
	// {
	// 	int count = 0;
	// 	for(int m=0; m<cp_exist.size(); m++){
	// 		auto tmp = m;
	// 		int i = tmp / dim0_offset;
	// 		tmp = tmp % dim0_offset;
	// 		int j = tmp / dim1_offset;
	// 		tmp = tmp % dim1_offset;
	// 		int k = tmp;
	// 		if(cp_exist[m]) count ++;
	// 	}
	// 	std::cout << count << std::endl;
	// 	std::cout << + cp_exist[4338498] << std::endl;
	// }
	std::cout << "start compression\n";
	for(int i=0; i<r1; i++){
		for(int j=0; j<r2; j++){
			for(int k=0; k<r3; k++){
				T required_eb = max_eb;
				// derive eb given 24 adjacent simplex
				for(int n=0; n<24; n++){
					bool in_mesh = true;
					for(int p=0; p<3; p++){
						// reversed order!
						if(!(in_range(i + index_offset[n][p][2], (int)r1) && in_range(j + index_offset[n][p][1], (int)r2) && in_range(k + index_offset[n][p][0], (int)r3))){
							in_mesh = false;
							break;
						}
					}
					if(in_mesh){
						int index = simplex_offset[n] + 6*(i*(r2-1)*(r3-1) + j*(r3-1) + k);
						// if(i*dim0_offset + j*dim1_offset + k == 987666){
						// 	std::cout << "simplex id = " << index << std::endl;
						// 	for(int i=0; i<3; i++){
						// 		std::cout << 987666 + offset[n][i] << " ";
						// 	}
						// 	std::cout << i*dim0_offset + j*dim1_offset + k << std::endl;
						// 	std::cout << "cp[" << index << "] = " << +cp_exist[index] << std::endl;
						// }
						if(cp_exist[index]){
							required_eb = 0;
							break;
						}
						required_eb = MINF(required_eb, derive_cp_abs_eb_sos_online<__int128, int64_t>(
							cur_U_pos[offset[n][0]], cur_U_pos[offset[n][1]], cur_U_pos[offset[n][2]], *cur_U_pos,
							cur_V_pos[offset[n][0]], cur_V_pos[offset[n][1]], cur_V_pos[offset[n][2]], *cur_V_pos,
							cur_W_pos[offset[n][0]], cur_W_pos[offset[n][1]], cur_W_pos[offset[n][2]], *cur_W_pos));
					}
				}
				// {
				// 	int id = i*dim0_offset + j*dim1_offset + k;
				// 	int target[4] = {725522, 987666, 988178, 988179};
				// 	if((id==target[0]) || (id==target[1]) || (id==target[2]) || (id==target[3])){
				// 		std::cout << id << ": eb = " << required_eb << std::endl;
				// 	}
				// }
				T abs_eb = required_eb;
				*eb_quant_index_pos = eb_exponential_quantize(abs_eb, base, log_of_base, threshold);
				// if(verbose) std::cout << std::endl << required_eb << " " << abs_eb << " unpred_size = " << unpred_data.size() << std::endl;
				if(abs_eb > 0){
					bool unpred_flag = false;
					T decompressed[3];
					// compress vector fields
					T * data_pos[3] = {cur_U_pos, cur_V_pos, cur_W_pos};
					for(int p=0; p<3; p++){
						T * cur_data_pos = data_pos[p];
						T cur_data = *cur_data_pos;
						// get adjacent data and perform Lorenzo
						/*
							d6	X
							d4	d5
							d2	d3
							d0	d1
						*/
						T d0 = (i && j && k) ? cur_data_pos[- dim0_offset - dim1_offset - 1] : 0;
						T d1 = (i && j) ? cur_data_pos[- dim0_offset - dim1_offset] : 0;
						T d2 = (i && k) ? cur_data_pos[- dim0_offset - 1] : 0;
						T d3 = (i) ? cur_data_pos[- dim0_offset] : 0;
						T d4 = (j && k) ? cur_data_pos[- dim1_offset - 1] : 0;
						T d5 = (j) ? cur_data_pos[- dim1_offset] : 0;
						T d6 = (k) ? cur_data_pos[- 1] : 0;
						T pred = d0 + d3 + d5 + d6 - d1 - d2 - d4;
						T diff = cur_data - pred;
						T quant_diff = std::abs(diff) / abs_eb + 1;
						if(quant_diff < capacity){
							quant_diff = (diff > 0) ? quant_diff : -quant_diff;
							int quant_index = (int)(quant_diff/2) + intv_radius;
							data_quant_index_pos[p] = quant_index;
							decompressed[p] = pred + 2 * (quant_index - intv_radius) * abs_eb; 
							// check original data
							if(std::abs(decompressed[p] - cur_data) >= required_eb){
								unpred_flag = true;
								break;
							}
						}
						else{
							unpred_flag = true;
							break;
						}
					}
					if(unpred_flag){
						*(eb_quant_index_pos ++) = 0;
						ptrdiff_t offset = cur_U_pos - U_fp;
						unpred_data.push_back(U[offset]);
						unpred_data.push_back(V[offset]);
						unpred_data.push_back(W[offset]);
					}
					else{
						eb_quant_index_pos ++;
						data_quant_index_pos += 3;
						*cur_U_pos = decompressed[0];
						*cur_V_pos = decompressed[1];
						*cur_W_pos = decompressed[2];
					}
				}
				else{
					// record as unpredictable data
					*(eb_quant_index_pos ++) = 0;
					ptrdiff_t offset = cur_U_pos - U_fp;
					unpred_data.push_back(U[offset]);
					unpred_data.push_back(V[offset]);
					unpred_data.push_back(W[offset]);
				}
				cur_U_pos ++, cur_V_pos ++, cur_W_pos ++;
			}
		}
	}
	free(U_fp);
	free(V_fp);
	free(W_fp);
	printf("offset eb_q, data_q, unpred: %ld %ld %ld\n", eb_quant_index_pos - eb_quant_index, data_quant_index_pos - data_quant_index, unpred_data.size());
	unsigned char * compressed = (unsigned char *) malloc(3*num_elements*sizeof(T));
	unsigned char * compressed_pos = compressed;
	write_variable_to_dst(compressed_pos, vector_field_scaling_factor);
	write_variable_to_dst(compressed_pos, base);
	write_variable_to_dst(compressed_pos, threshold);
	write_variable_to_dst(compressed_pos, intv_radius);
	size_t unpredictable_count = unpred_data.size();
	write_variable_to_dst(compressed_pos, unpredictable_count);
	write_array_to_dst(compressed_pos, (T_data *)&unpred_data[0], unpredictable_count);	
	size_t eb_quant_num = eb_quant_index_pos - eb_quant_index;
	write_variable_to_dst(compressed_pos, eb_quant_num);
	Huffman_encode_tree_and_data(2*1024, eb_quant_index, num_elements, compressed_pos);
	free(eb_quant_index);
	size_t data_quant_num = data_quant_index_pos - data_quant_index;
	write_variable_to_dst(compressed_pos, data_quant_num);
	Huffman_encode_tree_and_data(2*capacity, data_quant_index, data_quant_num, compressed_pos);
	printf("pos = %ld\n", compressed_pos - compressed);
	free(data_quant_index);
	compressed_size = compressed_pos - compressed;
	return compressed;	
}

template
unsigned char *
sz_compress_cp_preserve_sos_3d_online_fp(const float * U, const float * V, const float * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_sos_3d_online_fp(const double * U, const double * V, const double * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb);
