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

// default coordinates for tets in a cell
static const double default_coords[6][4][3] = {
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
  },
};

// compute offsets for simplex, index, and positions
static void 
compute_offset(ptrdiff_t dim0_offset, ptrdiff_t dim1_offset, ptrdiff_t cell_dim0_offset, ptrdiff_t cell_dim1_offset,
				int simplex_offset[24], int index_offset[24][3][3], int offset[24][3]){
	int * simplex_offset_pos = simplex_offset;
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
	for(int i=0; i<24; i++){
		for(int j=0; j<3; j++){
			for(int k=0; k<3; k++){
				index_offset[i][j][k] = coordinates[i][j][k] - coordinates[i][3][k];
			}
		}
	}
	for(int i=0; i<24; i++){
		for(int x=0; x<3; x++){
			offset[i][x] = (coordinates[i][x][0] - coordinates[i][3][0]) + (coordinates[i][x][1] - coordinates[i][3][1]) * dim1_offset + (coordinates[i][x][2] - coordinates[i][3][2]) * dim0_offset;
		}
	}	
}

template<typename T_fp>
static int 
check_cp(T_fp vf[4][3], int indices[4]){
	// robust critical point test
	bool succ = ftk::robust_critical_point_in_simplex3(vf, indices);
	if (!succ) return -1;
	return 1;
}

#define SINGULAR 0
#define STABLE_SOURCE 1
#define UNSTABLE_SOURCE 2
#define STABLE_REPELLING_SADDLE 3
#define UNSTABLE_REPELLING_SADDLE 4
#define STABLE_ATRACTTING_SADDLE  5
#define UNSTABLE_ATRACTTING_SADDLE  6
#define STABLE_SINK 7
#define UNSTABLE_SINK 8

template<typename T>
int
get_cp_type(const T X[4][3], const T U[4][3]){
	const T X_[3][3] = {
		{X[0][0] - X[3][0], X[1][0] - X[3][0], X[2][0] - X[3][0]}, 
		{X[0][1] - X[3][1], X[1][1] - X[3][1], X[2][1] - X[3][1]},
		{X[0][2] - X[3][2], X[1][2] - X[3][2], X[2][2] - X[3][2]}    
	};
	const T U_[3][3] = {
		{U[0][0] - U[3][0], U[1][0] - U[3][0], U[2][0] - U[3][0]}, 
		{U[0][1] - U[3][1], U[1][1] - U[3][1], U[2][1] - U[3][1]},
		{U[0][2] - U[3][2], U[1][2] - U[3][2], U[2][2] - U[3][2]}    
	};
	T inv_X_[3][3];
	ftk::matrix_inverse3x3(X_, inv_X_);
	T J[3][3];
	ftk::matrix3x3_matrix3x3_multiplication(inv_X_, U_, J);
	T P[4];
	ftk::characteristic_polynomial_3x3(J, P);
	std::complex<T> root[3];
	T disc = ftk::solve_cubic(P[2], P[1], P[0], root);
	if(fabs(disc) < std::numeric_limits<T>::epsilon()) return SINGULAR;
	int negative_real_parts = 0;
	for(int i=0; i<3; i++){
		negative_real_parts += (root[i].real() < 0);
	}
	switch(negative_real_parts){
		case 0:
			return (disc > 0) ? UNSTABLE_SOURCE : STABLE_SOURCE;
		case 1:
			return (disc > 0) ? UNSTABLE_REPELLING_SADDLE : STABLE_REPELLING_SADDLE;
		case 2:
			return (disc > 0) ? UNSTABLE_ATRACTTING_SADDLE : STABLE_ATRACTTING_SADDLE;
		case 3:
			return (disc > 0) ? UNSTABLE_SINK : STABLE_SINK;
		default:
			return SINGULAR;
	}
}

template<typename T_fp>
static int 
check_cp_type(const T_fp vf[4][3], const double v[4][3], const double X[4][3], const int indices[4]){
	// robust critical point test
	bool succ = ftk::robust_critical_point_in_simplex3(vf, indices);
	if (!succ) return -1;
	return get_cp_type(X, v);
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
	ptrdiff_t dim0_offset = r2*r3;
	ptrdiff_t dim1_offset = r3;
	ptrdiff_t cell_dim0_offset = (r2-1)*(r3-1);
	ptrdiff_t cell_dim1_offset = r3-1;
	int indices[4] = {0};
	// __int128 vf[4][3] = {0};
	int64_t vf[4][3] = {0};
	for(int i=0; i<r1-1; i++){
		for(int j=0; j<r2-1; j++){
			for(int k=0; k<r3-1; k++){
				bool verbose = false;
				// order (reserved, z->x):
				ptrdiff_t cell_offset = 6*(i*cell_dim0_offset + j*cell_dim1_offset + k);
				ptrdiff_t tmp_offset = 6*(i*dim0_offset + j*dim1_offset + k);
				// if(tmp_offset/6 == 4001374/6) verbose = true;
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
				// if(verbose){
				// 	auto offset = cell_offset + 3;
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
				// (ftk-3) 000, 010, 110, 111
				update_index(vf, indices, 1, i*dim0_offset + (j+1)*dim1_offset + k, U_fp, V_fp, W_fp);
				update_index(vf, indices, 2, i*dim0_offset + (j+1)*dim1_offset + k+1, U_fp, V_fp, W_fp);
				cp_exist[cell_offset + 4] = (check_cp(vf, indices) == 1);
				// if(verbose){
				// 	auto offset = cell_offset + 4;
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
				// (ftk-5) 000, 100, 110, 111
				update_index(vf, indices, 1, i*dim0_offset + j*dim1_offset + k+1, U_fp, V_fp, W_fp);
				cp_exist[cell_offset + 5] = (check_cp(vf, indices) == 1);
			}
		}
	}
	return cp_exist;	
}

template<typename T_data>
static inline void 
update_value(double v[4][3], int local_id, int global_id, const T_data * U, const T_data * V, const T_data * W){
	v[local_id][0] = U[global_id];
	v[local_id][1] = V[global_id];
	v[local_id][2] = W[global_id];
}

template<typename T_data, typename T_fp>
static vector<int> 
compute_cp_and_type(const T_fp * U_fp, const T_fp * V_fp, const T_fp * W_fp, const T_data * U, const T_data * V, const T_data * W, int r1, int r2, int r3){
	// check cp for all cells
	vector<int> cp_type(6*(r1-1)*(r2-1)*(r3-1), 0);
	ptrdiff_t dim0_offset = r2*r3;
	ptrdiff_t dim1_offset = r3;
	ptrdiff_t cell_dim0_offset = (r2-1)*(r3-1);
	ptrdiff_t cell_dim1_offset = r3-1;
	int indices[4] = {0};
	// __int128 vf[4][3] = {0};
	int64_t vf[4][3] = {0};
	double v[4][3] = {0};
	for(int i=0; i<r1-1; i++){
		for(int j=0; j<r2-1; j++){
			for(int k=0; k<r3-1; k++){
				bool verbose = false;
				// order (reserved, z->x):
				ptrdiff_t cell_offset = 6*(i*cell_dim0_offset + j*cell_dim1_offset + k);
				ptrdiff_t tmp_offset = 6*(i*dim0_offset + j*dim1_offset + k);
				// (ftk-0) 000, 001, 011, 111
				update_index(vf, indices, 0, i*dim0_offset + j*dim1_offset + k, U_fp, V_fp, W_fp);
				update_index(vf, indices, 1, (i+1)*dim0_offset + j*dim1_offset + k, U_fp, V_fp, W_fp);
				update_index(vf, indices, 2, (i+1)*dim0_offset + (j+1)*dim1_offset + k, U_fp, V_fp, W_fp);
				update_index(vf, indices, 3, (i+1)*dim0_offset + (j+1)*dim1_offset + (k+1), U_fp, V_fp, W_fp);
				for(int p=0; p<4; p++){
					v[p][0] = U[indices[p]];
					v[p][1] = V[indices[p]];
					v[p][2] = W[indices[p]];
				}
				cp_type[cell_offset] = check_cp_type(vf, v, default_coords[0], indices);
				// (ftk-2) 000, 010, 011, 111
				update_index(vf, indices, 1, i*dim0_offset + (j+1)*dim1_offset + k, U_fp, V_fp, W_fp);
				update_value(v, 1, indices[1], U, V, W); 
				cp_type[cell_offset + 1] = check_cp_type(vf, v, default_coords[1], indices);
				// (ftk-1) 000, 001, 101, 111
				update_index(vf, indices, 1, (i+1)*dim0_offset + j*dim1_offset + k, U_fp, V_fp, W_fp);
				update_value(v, 1, indices[1], U, V, W); 
				update_index(vf, indices, 2, (i+1)*dim0_offset + j*dim1_offset + k+1, U_fp, V_fp, W_fp);
				update_value(v, 2, indices[2], U, V, W); 
				cp_type[cell_offset + 2] = check_cp_type(vf, v, default_coords[2], indices);
				// (ftk-4) 000, 100, 101, 111
				update_index(vf, indices, 1, i*dim0_offset + j*dim1_offset + k+1, U_fp, V_fp, W_fp);
				update_value(v, 1, indices[1], U, V, W); 
				cp_type[cell_offset + 3] = check_cp_type(vf, v, default_coords[3], indices);
				// (ftk-3) 000, 010, 110, 111
				update_index(vf, indices, 1, i*dim0_offset + (j+1)*dim1_offset + k, U_fp, V_fp, W_fp);
				update_value(v, 1, indices[1], U, V, W); 
				update_index(vf, indices, 2, i*dim0_offset + (j+1)*dim1_offset + k+1, U_fp, V_fp, W_fp);
				update_value(v, 2, indices[2], U, V, W); 
				cp_type[cell_offset + 4] = check_cp_type(vf, v, default_coords[4], indices);
				// (ftk-5) 000, 100, 110, 111
				update_index(vf, indices, 1, i*dim0_offset + j*dim1_offset + k+1, U_fp, V_fp, W_fp);
				update_value(v, 1, indices[1], U, V, W); 
				cp_type[cell_offset + 5] = check_cp_type(vf, v, default_coords[5], indices);
			}
		}
	}
	return cp_type;	
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

std::ostream& operator<<(std::ostream& o, const __int128& x) {
    if (x == std::numeric_limits<__int128>::min()) return o << "-170141183460469231731687303715884105728";
    if (x < 0) return o << "-" << -x;
    if (x < 10) return o << (char)(x + '0');
    return o << x / 10 << (char)(x % 10 + '0');
}
/*
Tet x0, x1, x2, x3, derive cp-preserving eb for x3 given x0, x1, x2
using SoS method
*/
template<typename T>
static T 
derive_cp_abs_eb_sos_online(const T u0, const T u1, const T u2, const T u3, const T v0, const T v1, const T v2, const T v3, const T w0, const T w1, const T w2, const T w3, bool verbose=false){
	T M0 = - det_3_by_3(u1, u2, u3, v1, v2, v3, w1, w2, w3);
	T M1 = + det_3_by_3(u0, u2, u3, v0, v2, v3, w0, w2, w3);
	T M2 = - det_3_by_3(u0, u1, u3, v0, v1, v3, w0, w1, w3);
	T M3 = + det_3_by_3(u0, u1, u2, v0, v1, v2, w0, w1, w2);
	T M = M0 + M1 + M2 + M3;
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
	T denominator = abs(det_3_by_3(v0, v1, v2, w0, w1, w2, one, one, one)) + abs(det_3_by_3(u0, u1, u2, w0, w1, w2, one, one, one)) + abs(det_3_by_3(u0, u1, u2, v0, v1, v2, one, one, one)); 
	T eb = abs(M) / denominator;
	{
		// keep sign for replacing the three other vertices
		denominator = std::abs(det_2_by_2(v1, v2, w1, w2)) + std::abs(det_2_by_2(u1, u2, w1, w2)) + std::abs(det_2_by_2(u1, u2, v1, v2));
		if(denominator != 0){
			eb = MINF(eb, abs(M0) / denominator);
		}
		else return 0;
		denominator = std::abs(det_2_by_2(v0, v2, w0, w2)) + std::abs(det_2_by_2(u0, u2, w0, w2)) + std::abs(det_2_by_2(u0, u2, v0, v2));
		if(denominator != 0){
			eb = MINF(eb, abs(M1) / denominator);
		}
		else return 0;
		denominator = std::abs(det_2_by_2(v0, v1, w0, w1)) + std::abs(det_2_by_2(u0, u1, w0, w1)) + std::abs(det_2_by_2(u0, u1, v0, v1));
		if(denominator != 0){
			eb = MINF(eb, abs(M2) / denominator);
		}
		else return 0;
		// T cur_eb_0 = abs(M0)/(std::abs(det_2_by_2(v1, v2, w1, w2)) + std::abs(det_2_by_2(u1, u2, w1, w2)) + std::abs(det_2_by_2(u1, u2, v1, v2)));
		// T cur_eb_1 = abs(M1)/(std::abs(det_2_by_2(v0, v2, w0, w2)) + std::abs(det_2_by_2(u0, u2, w0, w2)) + std::abs(det_2_by_2(u0, u2, v0, v2)));
		// T cur_eb_2 = abs(M2)/(std::abs(det_2_by_2(v0, v1, w0, w1)) + std::abs(det_2_by_2(u0, u1, w0, w1)) + std::abs(det_2_by_2(u0, u1, v0, v1)));
		// eb = MINF(MINF(cur_eb_0, cur_eb_1), MINF(cur_eb_2, eb));
	}
	return eb;
}

template<typename T_acc, typename T>
static T 
derive_cp_abs_eb_sos_online_acc(const T u0, const T u1, const T u2, const T u3, const T v0, const T v1, const T v2, const T v3, const T w0, const T w1, const T w2, const T w3, bool verbose=false){
	T_acc M0 = - det_3_by_3<T_acc>(u1, u2, u3, v1, v2, v3, w1, w2, w3);
	// if(verbose){
	// 	std::cout << "M0 = " << M0 << std::endl;
	// }
	T_acc M1 = + det_3_by_3<T_acc>(u0, u2, u3, v0, v2, v3, w0, w2, w3);
	// if(verbose){
	// 	std::cout << "M1 = " << M1 << std::endl;
	// }
	T_acc M2 = - det_3_by_3<T_acc>(u0, u1, u3, v0, v1, v3, w0, w1, w3);
	// if(verbose){
	// 	std::cout << "M2 = " << M2 << std::endl;
	// }
	T_acc M3 = + det_3_by_3<T_acc>(u0, u1, u2, v0, v1, v2, w0, w1, w2);
	// if(verbose){
	// 	std::cout << "M3 = " << M3 << std::endl;
	// }
	T_acc M = M0 + M1 + M2 + M3;
	// if(verbose){
	// 	std::cout << u0 << " " << v0 << " " << w0 << std::endl;
	// 	std::cout << u1 << " " << v1 << " " << w1 << std::endl;
	// 	std::cout << u2 << " " << v2 << " " << w2 << std::endl;
	// 	std::cout << u3 << " " << v3 << " " << w3 << std::endl;
	// 	std::cout << "det = " << M << ": " << M0 << " " << M1 << " " << M2 << " " << M3 << std::endl;
	// }
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
	// if(verbose){
	// 	std::cout << "same_eb = " << same_eb << std::endl;
	// }
	if(same_eb != 0) return same_eb;
	// keep sign for the original simplex
	T one = 1;
	T_acc denominator = abs(det_3_by_3(v0, v1, v2, w0, w1, w2, one, one, one)) + abs(det_3_by_3(u0, u1, u2, w0, w1, w2, one, one, one)) + abs(det_3_by_3(u0, u1, u2, v0, v1, v2, one, one, one)); 
	T_acc eb = abs(M) / denominator;
	{
		// keep sign for replacing the three other vertices
		T_acc cur_eb_0 = abs(M0)/(std::abs(det_2_by_2(v1, v2, w1, w2)) + std::abs(det_2_by_2(u1, u2, w1, w2)) + std::abs(det_2_by_2(u1, u2, v1, v2)));
		T_acc cur_eb_1 = abs(M1)/(std::abs(det_2_by_2(v0, v2, w0, w2)) + std::abs(det_2_by_2(u0, u2, w0, w2)) + std::abs(det_2_by_2(u0, u2, v0, v2)));
		T_acc cur_eb_2 = abs(M2)/(std::abs(det_2_by_2(v0, v1, w0, w1)) + std::abs(det_2_by_2(u0, u1, w0, w1)) + std::abs(det_2_by_2(u0, u1, v0, v1)));
		// if(verbose){
		// 	T d1 = det_2_by_2(v1, v2, w1, w2);
		// 	T d2 = det_2_by_2(u1, u2, w1, w2);
		// 	T d3 = det_2_by_2(u1, u2, v1, v2);
		// 	std::cout << "denominator 1 = " << d1 << " " << d2 << " " << d3 << std::endl;
		// 	std::cout << "M0 = " << (T_acc) u3 * d1 - (T_acc) v3 * d2 + (T_acc) w3 * d3 << std::endl;
		// 	std::cout << "eb 1-3 = " << cur_eb_0 << " " << cur_eb_1 << " " << cur_eb_2 << std::endl;
		// }
		eb = MINF(MINF(cur_eb_0, cur_eb_1), MINF(cur_eb_2, eb));
	}
	return (T) MINF(eb, (T_acc) std::numeric_limits<int64_t>::max());
}

template<typename T, typename T_fp>
static int64_t 
convert_to_fixed_point(const T * U, const T * V, const T * W, size_t num_elements, T_fp * U_fp, T_fp * V_fp, T_fp * W_fp, T_fp& range, int type_bits=63){
	double vector_field_resolution = 0;
	int64_t vector_field_scaling_factor = 1;
	for (int i=0; i<num_elements; i++){
		double min_val = std::max(std::max(fabs(U[i]), fabs(V[i])), fabs(W[i]));
		vector_field_resolution = std::max(vector_field_resolution, min_val);
	}
	int vbits = std::ceil(std::log2(vector_field_resolution));
	int nbits = (type_bits - 5) / 3;
	vector_field_scaling_factor = 1 << (nbits - vbits);
	std::cerr << "resolution=" << vector_field_resolution 
	<< ", factor=" << vector_field_scaling_factor 
	<< ", nbits=" << nbits << ", vbits=" << vbits << ", shift_bits=" << nbits - vbits << std::endl;
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
	int simplex_offset[24];
	int index_offset[24][3][3];
	int offset[24][3];
	compute_offset(dim0_offset, dim1_offset, cell_dim0_offset, cell_dim1_offset, simplex_offset, index_offset, offset);
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
	// 	exit(0);
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
						if(cp_exist[index]){
							required_eb = 0;
							break;
						}
						// bool verbose = false;
						// int ftk_index = i*dim0_offset + j*dim1_offset + k;
						// if(index == 3987429){
						// 	std::cout << "\nsimplex " << n << std::endl;
						// 	std::cout << "cp_exist = " << +cp_exist[index] << std::endl;
						// 	std::cout << "actual cell index = " << index << ": ";
						// 	std::cout << i << " " << j << " " << k << std::endl;
						// 	verbose = true;
						// 	__int128 vf[4][3];
						// 	int indices[4];
						// 	int ftk_index = i*dim0_offset + j*dim1_offset + k;
						// 	for(int i=0; i<3; i++){
						// 		indices[i] = ftk_index + offset[n][i];
						// 	}
						// 	indices[3] = ftk_index;
						// 	for(int i=0; i<4; i++){
						// 		std::cout << indices[i] << " ";
						// 	}
						// 	std::cout << std::endl;
						// 	vf[0][0] = cur_U_pos[offset[n][0]], vf[1][0] = cur_U_pos[offset[n][1]], vf[2][0] = cur_U_pos[offset[n][2]], vf[3][0] = *cur_U_pos;
						// 	vf[0][1] = cur_V_pos[offset[n][0]], vf[1][1] = cur_V_pos[offset[n][1]], vf[2][1] = cur_V_pos[offset[n][2]], vf[3][1] = *cur_V_pos;
						// 	vf[0][2] = cur_W_pos[offset[n][0]], vf[1][2] = cur_W_pos[offset[n][1]], vf[2][2] = cur_W_pos[offset[n][2]], vf[3][2] = *cur_W_pos;
						// 	std::cout << "cp_exist before compression = " << check_cp(vf, indices) << std::endl;

						// }
						required_eb = MINF(required_eb, derive_cp_abs_eb_sos_online(
							cur_U_pos[offset[n][0]], cur_U_pos[offset[n][1]], cur_U_pos[offset[n][2]], *cur_U_pos,
							cur_V_pos[offset[n][0]], cur_V_pos[offset[n][1]], cur_V_pos[offset[n][2]], *cur_V_pos,
							cur_W_pos[offset[n][0]], cur_W_pos[offset[n][1]], cur_W_pos[offset[n][2]], *cur_W_pos));
						// required_eb = MINF(required_eb, derive_cp_abs_eb_sos_online_acc<__int128, int64_t>(
						// 	cur_U_pos[offset[n][0]], cur_U_pos[offset[n][1]], cur_U_pos[offset[n][2]], *cur_U_pos,
						// 	cur_V_pos[offset[n][0]], cur_V_pos[offset[n][1]], cur_V_pos[offset[n][2]], *cur_V_pos,
						// 	cur_W_pos[offset[n][0]], cur_W_pos[offset[n][1]], cur_W_pos[offset[n][2]], *cur_W_pos));
						// if(verbose) exit(0);
					}
				}
				// {
				// 	int id = i*dim0_offset + j*dim1_offset + k;
				// 	int target[4] = {666895, 666896, 929040, 929552};
				// 	if((id==target[0]) || (id==target[1]) || (id==target[2]) || (id==target[3])){
				// 		std::cout << id << ": eb = " << required_eb << std::endl;
				// 	}
				// }
				T abs_eb = required_eb;
				*eb_quant_index_pos = eb_exponential_quantize(abs_eb, base, log_of_base, threshold);
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
	Huffman_encode_tree_and_data(2*1024, eb_quant_index, eb_quant_num, compressed_pos);
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

template<typename T_data>
unsigned char *
sz_compress_cp_preserve_sos_3d_online_fp_spec_eb(const T_data * U, const T_data * V, const T_data * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb){
	std::cout << "sz_compress_cp_preserve_sos_3d_online_fp_spec_eb" << std::endl;
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
	int simplex_offset[24];
	int index_offset[24][3][3];
	int offset[24][3];
	compute_offset(dim0_offset, dim1_offset, cell_dim0_offset, cell_dim1_offset, simplex_offset, index_offset, offset);
	T * cur_U_pos = U_fp;
	T * cur_V_pos = V_fp;
	T * cur_W_pos = W_fp;
	T threshold = 1;
	// check cp for all cells
	std::cout << "start cp checking\n";
	vector<bool> cp_exist = compute_cp(U_fp, V_fp, W_fp, r1, r2, r3);
	std::cout << "start compression\n";
	// intermediate variable
	T decompressed[3];
	T data_ori[3];
	T pred[3];
	T pred_residue[3];
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
						if(cp_exist[index]){
							required_eb = 0;
							break;
						}
						required_eb = MINF(required_eb, derive_cp_abs_eb_sos_online(
							cur_U_pos[offset[n][0]], cur_U_pos[offset[n][1]], cur_U_pos[offset[n][2]], *cur_U_pos,
							cur_V_pos[offset[n][0]], cur_V_pos[offset[n][1]], cur_V_pos[offset[n][2]], *cur_V_pos,
							cur_W_pos[offset[n][0]], cur_W_pos[offset[n][1]], cur_W_pos[offset[n][2]], *cur_W_pos));
						// required_eb = MINF(required_eb, derive_cp_abs_eb_sos_online_acc<__int128, int64_t>(
						// 	cur_U_pos[offset[n][0]], cur_U_pos[offset[n][1]], cur_U_pos[offset[n][2]], *cur_U_pos,
						// 	cur_V_pos[offset[n][0]], cur_V_pos[offset[n][1]], cur_V_pos[offset[n][2]], *cur_V_pos,
						// 	cur_W_pos[offset[n][0]], cur_W_pos[offset[n][1]], cur_W_pos[offset[n][2]], *cur_W_pos));
						// if(verbose) exit(0);
					}
				}
				T abs_eb = required_eb;
				// relax error bound
				abs_eb = relax_eb(abs_eb, (T) 8);
				abs_eb = MINF(abs_eb, max_eb);
				*eb_quant_index_pos = eb_exponential_quantize(abs_eb, base, log_of_base, threshold);
				if(abs_eb > 0){
					bool unpred_flag = false;
					// compress vector fields
					T * data_pos[3] = {cur_U_pos, cur_V_pos, cur_W_pos};
					for(int p=0; p<3; p++){
						T * cur_data_pos = data_pos[p];
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
						data_ori[p] = *cur_data_pos;
						pred[p] = d0 + d3 + d5 + d6 - d1 - d2 - d4;
						pred_residue[p] = data_ori[p] - pred[p];
					}
					bool repeated = false;
					while(true){
						bool verification_flag = true;
						// quantize using relaxed eb
						for(int p=0; p<3; p++){
							T diff = pred_residue[p];
							T quant_diff = std::abs(diff) / abs_eb + 1;
							if(quant_diff < capacity){
								quant_diff = (diff > 0) ? quant_diff : -quant_diff;
								int quant_index = (int)(quant_diff/2) + intv_radius;
								data_quant_index_pos[p] = quant_index;
								decompressed[p] = pred[p] + 2 * (quant_index - intv_radius) * abs_eb; 
								// check original data
								if(std::abs(decompressed[p] - data_ori[p]) >= required_eb){
									verification_flag = false;
									break;
								}
							}
							else{
								unpred_flag = true;
								verification_flag = true;
								break;
							}					
						}			
						if(verification_flag) break;
						if(repeated){
							// special case when decompressed[p] - data_ori[p] == required_eb
							// will lead to finite loop
							unpred_flag = true;
							break;
						}
						abs_eb = restrict_eb(abs_eb);
						if(abs_eb < required_eb){
							abs_eb = required_eb;
							repeated = true;
						}
						*eb_quant_index_pos = eb_exponential_quantize(abs_eb, base, log_of_base, threshold);
						if(abs_eb == 0){
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
	Huffman_encode_tree_and_data(2*1024, eb_quant_index, eb_quant_num, compressed_pos);
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
sz_compress_cp_preserve_sos_3d_online_fp_spec_eb(const float * U, const float * V, const float * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_sos_3d_online_fp_spec_eb(const double * U, const double * V, const double * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb);

template<typename T_data>
unsigned char *
sz_compress_cp_preserve_sos_3d_online_fp_spec_exec_fn(const T_data * U, const T_data * V, const T_data * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb, double max_factor){
	std::cout << "sz_compress_cp_preserve_sos_3d_online_fp_spec_exec_fn" << std::endl;
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
	int simplex_offset[24];
	int index_offset[24][3][3];
	int offset[24][3];
	compute_offset(dim0_offset, dim1_offset, cell_dim0_offset, cell_dim1_offset, simplex_offset, index_offset, offset);
	T * cur_U_pos = U_fp;
	T * cur_V_pos = V_fp;
	T * cur_W_pos = W_fp;
	T threshold = 1;
	// check cp for all cells
	std::cout << "start cp checking\n";
	vector<bool> cp_exist = compute_cp(U_fp, V_fp, W_fp, r1, r2, r3);
	std::cout << "start compression\n";
	// intermediate variable
	T decompressed[3];
	int indices[4];
	T vf[4][3];
	for(int i=0; i<r1; i++){
		for(int j=0; j<r2; j++){
			for(int k=0; k<r3; k++){
				bool unpred_flag = false;
				bool verification_flag = false;
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
						if(cp_exist[index]){
							unpred_flag = true;
							verification_flag = true;
							break;
						}
					}
				}
				T abs_eb = max_eb;
				// compress data and then verify
				while(!verification_flag){
					*eb_quant_index_pos = eb_exponential_quantize(abs_eb, base, log_of_base, threshold);
					unpred_flag = false;
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
						}
						else{
							unpred_flag = true;
							break;
						}
					}
					if(unpred_flag) break;
					// verify cp in 24 adjacent triangles
					verification_flag = true;
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
							int vertex_index = i*dim0_offset + j*dim1_offset + k;
							for(int p=0; p<3; p++){
								indices[p] = vertex_index + offset[n][p];
							}
							indices[3] = vertex_index;
							// TODO: change stragegy and consider types
							for(int p=0; p<3; p++){
								vf[p][0] = U_fp[indices[p]];
								vf[p][1] = V_fp[indices[p]];
								vf[p][2] = W_fp[indices[p]];
							}
							vf[3][0] = decompressed[0], vf[3][1] = decompressed[1], vf[3][2] = decompressed[2];
							bool decompressed_has_cp = (check_cp(vf, indices) == 1);
							if(decompressed_has_cp){
								verification_flag = false;
								break;
							}
						}
					}
					// relax error bound
					abs_eb = restrict_eb(abs_eb);
					if((!verification_flag) && (abs_eb <= max_eb * 1.0/max_factor)){
						unpred_flag = true;
						verification_flag = true;					
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
	Huffman_encode_tree_and_data(2*1024, eb_quant_index, eb_quant_num, compressed_pos);
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
sz_compress_cp_preserve_sos_3d_online_fp_spec_exec_fn(const float * U, const float * V, const float * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb, double max_factor);

template
unsigned char *
sz_compress_cp_preserve_sos_3d_online_fp_spec_exec_fn(const double * U, const double * V, const double * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb, double max_factor);

template<typename T_data, typename T_fp>
static inline T_data convert_fp_to_float(T_fp fp, T_fp vector_field_scaling_factor){
	return fp * (T_data) 1.0 / vector_field_scaling_factor;
}
template<typename T_data>
unsigned char *
sz_compress_cp_preserve_sos_3d_online_fp_spec_exec_all(const T_data * U, const T_data * V, const T_data * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb, double max_factor){
	std::cout << "sz_compress_cp_preserve_sos_3d_online_fp_spec_exec_all" << std::endl;
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
	int simplex_offset[24];
	int index_offset[24][3][3];
	int offset[24][3];
	compute_offset(dim0_offset, dim1_offset, cell_dim0_offset, cell_dim1_offset, simplex_offset, index_offset, offset);
	double coordinates_d[24][4][3];
	for(int i=0; i<24; i++){
		for(int j=0; j<4; j++){
			for(int k=0; k<3; k++){
				coordinates_d[i][j][k] = coordinates[i][j][k];
			}
		}
	}
	T * cur_U_pos = U_fp;
	T * cur_V_pos = V_fp;
	T * cur_W_pos = W_fp;
	// dec_data
	T_data * dec_U = (T_data *) malloc(num_elements*sizeof(T_data));
	T_data * dec_V = (T_data *) malloc(num_elements*sizeof(T_data));
	T_data * dec_W = (T_data *) malloc(num_elements*sizeof(T_data));
	memcpy(dec_U, U, num_elements*sizeof(T_data));
	memcpy(dec_V, V, num_elements*sizeof(T_data));
	memcpy(dec_W, W, num_elements*sizeof(T_data));
	T threshold = 1;
	// check cp for all cells
	std::cout << "start cp checking\n";
	vector<int> cp_type = compute_cp_and_type(U_fp, V_fp, W_fp, U, V, W, r1, r2, r3);
	std::cout << "start compression\n";
	for(int i=0; i<r1; i++){
		for(int j=0; j<r2; j++){
			for(int k=0; k<r3; k++){
				T abs_eb = max_eb;
				bool unpred_flag = false;
				bool verification_flag = false;
				T decompressed[3];
				// compress data and then verify
				while(!verification_flag){
					*eb_quant_index_pos = eb_exponential_quantize(abs_eb, base, log_of_base, threshold);
					unpred_flag = false;
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
						}
						else{
							unpred_flag = true;
							break;
						}
					}
					if(unpred_flag) break;
					// verify cp in six adjacent triangles
					verification_flag = true;
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
							int indices[4];
							T vf[4][3];
							int vertex_index = i*dim0_offset + j*dim1_offset + k;
							for(int p=0; p<3; p++){
								indices[p] = vertex_index + offset[n][p];
							}
							indices[3] = vertex_index;
							for(int p=0; p<3; p++){
								vf[p][0] = U_fp[indices[p]];
								vf[p][1] = V_fp[indices[p]];
								vf[p][2] = W_fp[indices[p]];
							}
							vf[3][0] = decompressed[0], vf[3][1] = decompressed[1], vf[3][2] = decompressed[2];
							double v[4][3];
							// use decompressed/original data for other vertices
							for(int p=0; p<3; p++){
								v[p][0] = dec_U[indices[p]];
								v[p][1] = dec_V[indices[p]];
								v[p][2] = dec_W[indices[p]];
							}
							// compute decompressed data for current vertex
							for(int p=0; p<3; p++){
								v[3][p] = convert_fp_to_float<T_data>(decompressed[p], vector_field_scaling_factor);
							}
							int decompressed_cp_type = check_cp_type(vf, v, coordinates_d[n], indices);
							int cell_index = simplex_offset[n] + 6*(i*cell_dim0_offset + j*cell_dim1_offset + k);
							if(decompressed_cp_type != cp_type[cell_index]){
								verification_flag = false;
								break;
							}
						}
					}
					// relax error bound
					abs_eb = restrict_eb(abs_eb);
					if((!verification_flag) && (abs_eb <= max_eb * 1.0/max_factor)){
						unpred_flag = true;
						verification_flag = true;					
					}
				}
				ptrdiff_t offset = cur_U_pos - U_fp;
				if(unpred_flag){
					*(eb_quant_index_pos ++) = 0;
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
					dec_U[offset] = convert_fp_to_float<T_data>(decompressed[0], vector_field_scaling_factor);
					dec_V[offset] = convert_fp_to_float<T_data>(decompressed[1], vector_field_scaling_factor);
					dec_W[offset] = convert_fp_to_float<T_data>(decompressed[2], vector_field_scaling_factor);
				}
				cur_U_pos ++, cur_V_pos ++, cur_W_pos ++;
			}
		}
	}
	free(dec_U);
	free(dec_V);
	free(dec_W);
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
	Huffman_encode_tree_and_data(2*1024, eb_quant_index, eb_quant_num, compressed_pos);
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
sz_compress_cp_preserve_sos_3d_online_fp_spec_exec_all(const float * U, const float * V, const float * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb, double max_factor);

template
unsigned char *
sz_compress_cp_preserve_sos_3d_online_fp_spec_exec_all(const double * U, const double * V, const double * W, size_t r1, size_t r2, size_t r3, size_t& compressed_size, bool transpose, double max_pwr_eb, double max_factor);
