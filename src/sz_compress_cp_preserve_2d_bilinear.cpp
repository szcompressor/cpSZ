#include "sz_cp_preserve_utils.hpp"
#include "sz_compress_3d.hpp"
#include "sz_compress_cp_preserve_2d.hpp"
#include "sz_def.hpp"
#include "sz_compression_utils.hpp"

// maximal error bound to keep the sign of A*(1 + e_1) + B*(1 + e_2) + C
template<typename T>
static inline double max_eb_to_keep_sign_2d_online(const T A, const T B, const T C=0){
	double fabs_sum = (fabs(A) + fabs(B));
	if(fabs_sum == 0) return 0;
	return fabs(A + B + C) / fabs_sum;
}

// maximal error bound to keep the sign of 
// a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f
// when f < 0, c = +/- 2ab
static inline double max_eb_to_keep_sign_2d_online_lt0(const double a, const double b, const double c, const double d, const double e, const double f){
	double eb = 1;
	// check four corners
	{
		// [-1, -1]
		double tmp_a = a + b + c;
		double tmp_b = - d - e;
		double tmp_c = f;
		eb = MIN(eb, (-tmp_b + sqrt(tmp_b*tmp_b - 4*tmp_a*tmp_c))/(2*tmp_a));
	}
	{
		// [-1, 1]
		double tmp_a = a + b - c;
		double tmp_b = - d + e;
		double tmp_c = f;
		eb = MIN(eb, (-tmp_b + sqrt(tmp_b*tmp_b - 4*tmp_a*tmp_c))/(2*tmp_a));
	}
	{
		// [1, -1]
		double tmp_a = a + b - c;
		double tmp_b = d - e;
		double tmp_c = f;
		eb = MIN(eb, (-tmp_b + sqrt(tmp_b*tmp_b - 4*tmp_a*tmp_c))/(2*tmp_a));
	}
	{
		// [1, 1]
		double tmp_a = a + b + c;
		double tmp_b = d + e;
		double tmp_c = f;
		eb = MIN(eb, (-tmp_b + sqrt(tmp_b*tmp_b - 4*tmp_a*tmp_c))/(2*tmp_a));
	}			
	return eb;
}

// maximal error bound to keep the sign of 
// a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f
// when f > 0, c = +/- 2ab
// = a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f
// > 0 - (|d| + |e|) e1 + f
// > 0
static inline double max_eb_to_keep_sign_2d_online_gt0(const double a, const double b, const double c, const double d, const double e, const double f){
	return f / (fabs(d) + fabs(e));
}
/*
x1 - x2    rotate    x2 - x3
|    |     ------->  |    |
x0 - x3              x1 - x0
bilinear cell x0, x1, x2, x3, derive cp-preserving eb for x3 given x0, x1, x2
rotate to make (u3, v3) appear only once in A, B, C, D
*/
static double 
derive_cp_eb_bilinear_online(const double u0, const double u1, const double u2, const double u3, const double v0, const double v1, const double v2, const double v3, bool verbose=false){
	// f(u) = A0 xy + B0 x + C0 y + D0
	// f(v) = A1 xy + B1 x + C1 y + D1
	// solve 
	double A0 = u3 - u0 - u2 + u1;
	double B0 = u0 - u1;
	double C0 = u2 - u1;
	double D0 = u1;

	double A1 = v3 - v0 - v2 + v1;
	double B1 = v0 - v1;
	double C1 = v2 - v1;
	double D1 = v1;

	// M[x, 1]^T = lambda [x, 1]^T
	// 	1/(A0C1 - A1C0)	[-B0C1 + B1C0, -C1D0 + C0D1][x]	= y[x]
	// 									[ A1B0 - A0B1,  A1D0 - A0D1][1]		 [1]
	double A0C1_minus_A1C0 = A0*C1 - A1*C0;
	double M[4] = {- B0*C1 + B1*C0, - C1*D0 + C0*D1, A1*B0 - A0*B1,  A1*D0 - A0*D1};
	for(int i=0; i<4; i++){
		M[i] /= A0C1_minus_A1C0;
	}
	// original determinant
	double dM = M[0] * M[3] - M[1] * M[2];
	// original trace
	double tM = M[0] + M[3];

	// A0C1 - A1C0 = - u1v0 + u2v0 + u0v1 - u3v1 - u0v2 + u3v2 + u1v3 - u2v3

	// note: not dividing A0C1 - A1C0 in M
	// M[0, 0] = - u1v0 + u2v0 + u0v1 - u2v1 - u0v2 + u1v2
	// M[0, 1] = u2v1 - u1v2
	// M[1, 0] = u2v0 - u3v0 - u2v1 + u3v1 - u0v2 + u1v2 + u0v3 - u1v3
	// M[1, 1] = -u1v0 + u0v1 + u2v1 -u3v1 

	// tr(M) = -2 u1v0 + u2v0 + 2 u0v1 - u3v1 - u0v2 + u1v3
  // det(M) = u1u1v0v0 - u1u2v0v0 - 2 u0u1v0v1 + u0u2v0v1 + u1u3v0v1 + u0u0v1v1 
  //					- u0u3v1v1 + u0u1v0v2 - u1u3v0v2 - u0u0v1v2 + u0u3v1v2 - u1u1v0v3 
  //						+ u1u2v0v3 + u0u1v1v3 - u0u2v1v3
	double eb = 1;
	if(tM * tM - 4*dM == 0){
		// 1 double root
		return 0;
	}
	else if(tM * tM - 4*dM < 0) {
		// keep sign of delta
		// keep sign of trace_M * trace_M - 4 * determinant_M
		// tr(M)^2 - 4 det(M) = 
		// 	u2u2v0v0 - 2 u2u3v0v1 + u3u3v1v1 - 2 u0u2v0v2 + 4 u1u3v0v2 - 2 u0u3v1v2
		// +	u0u0v2v2 - 2 u1u2v0v3 + 4 u0u2v1v3 - 2 u1u3v1v3 - 2 u0u1v2v3 + u1u1v3v3
		// = a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f = 0
		double a = u3*u3*v1*v1; // (1+e1)u3(1+e1)u3v1v1
		double b = u1*u1*v3*v3; // u1u1(1+e2)v3(1+e2)v3
		double c = - 2*u1*u3*v1*v3; // -2 u1u3(1+e1)v1v3(1+e2)
		double d = - 2*u2*u3*v0*v1 + 2*u3*u3*v1*v1 + 4*u1*u3*v0*v2 - 2*u0*u3*v1*v2 - 2*u1*u3*v1*v3; // - 2 u2u3v0v1 + (1+e1)u3(1+e1)u3v1v1 + 4 u1u3v0v2 - 2 u0u3v1v2 - 2 u1(1+e1)u3v1(1+e2)v3 
		double e = - 2*u1*u2*v0*v3 + 4*u0*u2*v1*v3 - 2*u1*u3*v1*v3 - 2*u0*u1*v2*v3 + 2*u1*u1*v3*v3; // - 2 u1u2v0v3 + 4 u0u2v1v3 - 2 u1u3v1v3 - 2 u0u1v2v3 + u1u1v3v3
		double f = u2*u2*v0*v0 - 2*u2*u3*v0*v1 + u3*u3*v1*v1 - 2*u0*u2*v0*v2 + 4*u1*u3*v0*v2 - 2*u0*u3*v1*v2
		 +	u0*u0*v2*v2 - 2*u1*u2*v0*v3 + 4*u0*u2*v1*v3 - 2*u1*u3*v1*v3 - 2*u0*u1*v2*v3 + u1*u1*v3*v3;
		// keep sign of A0C1 - A1C0 = - u1v0 + u2v0 + u0v1 - u3v1 - u0v2 + u3v2 + u1v3 - u2v3
		{
			eb = MIN(eb, max_eb_to_keep_sign_2d_online(-u3*v1 + u3*v2, u1*v3 - u2*v3, -u1*v0 + u2*v0 + u0*v1 - u0*v2));
		}
		eb = MIN(eb, max_eb_to_keep_sign_2d_online_lt0(a, b, c, d, e, f));
	}		
	else{
		// return 0;
		// tM * tM - 4*dM > 0
		// keep sign of delta
		{
			double a = u3*u3*v1*v1; // (1+e1)u3(1+e1)u3v1v1
			double b = u1*u1*v3*v3; // u1u1(1+e2)v3(1+e2)v3
			double c = - 2*u1*u3*v1*v3; // -2 u1u3(1+e1)v1v3(1+e2)
			double d = - 2*u2*u3*v0*v1 + 2*u3*u3*v1*v1 + 4*u1*u3*v0*v2 - 2*u0*u3*v1*v2 - 2*u1*u3*v1*v3; // - 2 u2u3v0v1 + (1+e1)u3(1+e1)u3v1v1 + 4 u1u3v0v2 - 2 u0u3v1v2 - 2 u1(1+e1)u3v1(1+e2)v3 
			double e = - 2*u1*u2*v0*v3 + 4*u0*u2*v1*v3 - 2*u1*u3*v1*v3 - 2*u0*u1*v2*v3 + 2*u1*u1*v3*v3; // - 2 u1u2v0v3 + 4 u0u2v1v3 - 2 u1u3v1v3 - 2 u0u1v2v3 + u1u1v3v3
			double f = u2*u2*v0*v0 - 2*u2*u3*v0*v1 + u3*u3*v1*v1 - 2*u0*u2*v0*v2 + 4*u1*u3*v0*v2 - 2*u0*u3*v1*v2
			 +	u0*u0*v2*v2 - 2*u1*u2*v0*v3 + 4*u0*u2*v1*v3 - 2*u1*u3*v1*v3 - 2*u0*u1*v2*v3 + u1*u1*v3*v3;
			// keep sign of A0C1 - A1C0 = - u1v0 + u2v0 + u0v1 - u3v1 - u0v2 + u3v2 + u1v3 - u2v3
			{
				eb = MIN(eb, max_eb_to_keep_sign_2d_online(-u3*v1 + u3*v2, u1*v3 - u2*v3, -u1*v0 + u2*v0 + u0*v1 - u0*v2));
			}
			// = a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f
			// > 0 - (|d| + |e|) e1 + f
			// > 0
			eb = MIN(eb, max_eb_to_keep_sign_2d_online_gt0(a, b, c, d, e, f));
		}
		if((dM > 0) && (1 - tM + dM > 0)){
			// include 2 cases
			// 1) have roots but no root in [0, 1] case 1
			// f(0) > 0, f(1) > 0, Tr(M)/2 not in [0, 1]
			// 2) two roots in [0, 1]
			// f(0) > 0, f(1) > 0, Tr(M)/2 in [0, 1]
		  // keep sign of dM
		  {
			  double a = u1*u3*v0*v1 - u0*u3*v1*v1 - u1*u3*v0*v2 + u0*u3*v1*v2;
			  double b = - u1*u1*v0*v3 + u1*u2*v0*v3 + u0*u1*v1*v3 - u0*u2*v1*v3;
			  double c = u1*u1*v0*v0 - u1*u2*v0*v0 - 2*u0*u1*v0*v1 + u0*u2*v0*v1 + u0*u0*v1*v1 + u0*u1*v0*v2 - u0*u0*v1*v2; 
			  eb = MIN(eb, max_eb_to_keep_sign_2d_online(a, b, c));		  	
		  }
		  // tr(M) = (-2 u1v0 + u2v0 + 2 u0v1 - u3v1 - u0v2 + u1v3) / (- u1v0 + u2v0 + u0v1 - u3v1 - u0v2 + u3v2 + u1v3 - u2v3)
		  {
			  // keep sign of tM 
			  double a = - u3*v1;
			  double b = u1*v3;
			  double c = -2*u1*v0 + u2*v0 + 2*u0*v1 - u0*v2;
			  eb = MIN(eb, max_eb_to_keep_sign_2d_online(a, b, c));
			  // keep sign of tr(M) - 2(A0C1 - A1C0)
			  a -= 2*(- u3*v1 + u3*v2);
			  b -= 2*(u1*v3 - u2*v3);
			  c -= 2*(- u1*v0 + u2*v0 + u0*v1 - u0*v2);
			  eb = MIN(eb, max_eb_to_keep_sign_2d_online(a, b, c));
		  }
		  // keep sign of 1 - tM + dM
		  // (-2 u1v0 + u2v0 + 2 u0v1 - u3v1 - u0v2 + u1v3) * (- u1v0 + u2v0 + u0v1 - u3v1 - u0v2 + u3v2 + u1v3 - u2v3)
		  // = ?
		  // (- u1v0 + u2v0 + u0v1 - u3v1 - u0v2 + u3v2 + u1v3 - u2v3) ** 2
		  // = ?
		  {
		  	// double a = ;
		  	// double b = ;
		  	// double c = ;
		  	// double d = ;
		  	// double e = ;
		  	// double f = ;
		  	// if(f == 0) return 0;
		  	// else if(f > 0){
					// eb = MIN(eb, f / (fabs(d) + fabs(e)));
		  	// }
		  	// else{

		  	// }
		  }
		}
		else if((dM < 0) && (1 - tM + dM < 0)){
			return 0;
			// have roots but no root in [0, 1] case 2
			// f(0) < 0, f(1) < 0
			// keep sign of dM
		  double a = u1*u3*v0*v1 - u0*u3*v1*v1 - u1*u3*v0*v2 + u0*u3*v1*v2;
		  double b = - u1*u1*v0*v3 + u1*u2*v0*v3 + u0*u1*v1*v3 - u0*u2*v1*v3;
		  double c = u1*u1*v0*v0 - u1*u2*v0*v0 - 2*u0*u1*v0*v1 + u0*u2*v0*v1 + u0*u0*v1*v1 + u0*u1*v0*v2 - u0*u0*v1*v2; 
		  eb = MIN(eb, max_eb_to_keep_sign_2d_online(a, b, c));
		  // keep sign of 1 - tM + dM
		  a -= - u3*v1;
		  b -= u1*v3;
		  c -= -2*u1*v0 + u2*v0 + 2*u0*v1 - u0*v2;
		  c += 1;
		  eb = MIN(eb, max_eb_to_keep_sign_2d_online(a, b, c));
		}
		else{
			return 0;
			// dM * (1 - tM + dM) < 0
			// one root in [0, 1]
			// keep sign of dM
		  double a = u1*u3*v0*v1 - u0*u3*v1*v1 - u1*u3*v0*v2 + u0*u3*v1*v2;
		  double b = - u1*u1*v0*v3 + u1*u2*v0*v3 + u0*u1*v1*v3 - u0*u2*v1*v3;
		  double c = u1*u1*v0*v0 - u1*u2*v0*v0 - 2*u0*u1*v0*v1 + u0*u2*v0*v1 + u0*u0*v1*v1 + u0*u1*v0*v2 - u0*u0*v1*v2; 
		  eb = MIN(eb, max_eb_to_keep_sign_2d_online(a, b, c));
		  // keep sign of 1 - tM + dM
		  // tr(M) = -2 u1v0 + u2v0 + 2 u0v1 - u3v1 - u0v2 + u1v3
		  a -= - u3*v1;
		  b -= u1*v3;
		  c -= -2*u1*v0 + u2*v0 + 2*u0*v1 - u0*v2;
		  c += 1;
		  eb = MIN(eb, max_eb_to_keep_sign_2d_online(a, b, c));
		}
	}
	return eb;
}

template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_bilinear_online_log(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){

	size_t num_elements = r1 * r2;
	size_t sign_map_size = (num_elements - 1)/8 + 1;
	unsigned char * sign_map_compressed = (unsigned char *) malloc(2*sign_map_size);
	unsigned char * sign_map_compressed_pos = sign_map_compressed;
	unsigned char * sign_map = (unsigned char *) malloc(num_elements*sizeof(unsigned char));
	// Note the convert function has address auto increment
	T * log_U = log_transform(U, sign_map, num_elements);
	convertIntArray2ByteArray_fast_1b_to_result_sz(sign_map, num_elements, sign_map_compressed_pos);
	T * log_V = log_transform(V, sign_map, num_elements);
	convertIntArray2ByteArray_fast_1b_to_result_sz(sign_map, num_elements, sign_map_compressed_pos);
	free(sign_map);

	T * decompressed_U = (T *) malloc(num_elements*sizeof(T));
	memcpy(decompressed_U, U, num_elements*sizeof(T));
	T * decompressed_V = (T *) malloc(num_elements*sizeof(T));
	memcpy(decompressed_V, V, num_elements*sizeof(T));

	int * eb_quant_index = (int *) malloc(num_elements*sizeof(int));
	int * data_quant_index = (int *) malloc(2*num_elements*sizeof(int));
	int * eb_quant_index_pos = eb_quant_index;
	int * data_quant_index_pos = data_quant_index;
	// next, row by row
	const int base = 2;
	const double log_of_base = log2(base);
	const int capacity = 65536;
	const int intv_radius = (capacity >> 1);
	unpred_vec<T> unpred_data;
	// offsets to get six adjacent square mesh indices
	// T8 rolls back to T0
	/*
	|	T3	T4	T5
	y	T2	X 	T6
	|	T1	T0  T7
		-	x 	-
	*/
	const int offsets[9] = {
		-(int)r2, -(int)r2 - 1, -1, (int)r2 - 1, (int)r2, (int)r2+1, 1, -(int)r2 + 1, -(int)r2
	};
	// relative indices for the four squares
	// we always put the current index as (x3, y3)
	const T x[4][4] = {
		{1, 0, 0, 1}, // T0_T1_T2_X
		{0, 0, 1, 1}, // T2_T3_T4_X
		{0, 1, 1, 0}, // T4_T5_T6_X
		{1, 1, 0, 0}, // T6_T7_T0_X
	};
	const T y[4][4] = {
		{0, 0, 1, 1},
		{0, 1, 1, 0},
		{1, 1, 0, 0},
		{1, 0, 0, 1},
	};
  // compute offset of index for the surrounding cells	
	int index_offset[4][3][2];
	for(int i=0; i<4; i++){
		for(int j=0; j<3; j++){
			index_offset[i][j][0] = x[i][j] - x[i][3];
			index_offset[i][j][1] = y[i][j] - y[i][3];
		}
	}
	T * cur_log_U_pos = log_U;
	T * cur_log_V_pos = log_V;
	T * cur_U_pos = decompressed_U;
	T * cur_V_pos = decompressed_V;
	for(int i=0; i<r1; i++){
		// printf("start %d row\n", i);
		for(int j=0; j<r2; j++){
			double required_eb = max_pwr_eb;
			// derive eb given four adjacent cells
			for(int k=0; k<4; k++){
				bool in_mesh = true;
				for(int p=0; p<3; p++){
					// reserved order!
					// note: x and j are on r2, y and i are on r1
					if(!(in_range(i + index_offset[k][p][1], (int)r1) && in_range(j + index_offset[k][p][0], (int)r2))){
						in_mesh = false;
						break;
					}
				}
				if(in_mesh){
					required_eb = MIN(required_eb, derive_cp_eb_bilinear_online(cur_U_pos[offsets[2*k]], cur_U_pos[offsets[2*k+1]], cur_U_pos[offsets[2*k+2]], cur_U_pos[0],
						cur_V_pos[offsets[2*k]], cur_V_pos[offsets[2*k+1]], cur_V_pos[offsets[2*k+2]], cur_V_pos[0]));
				}
			}
			// eb[index_eb++] = required_eb;
			if((required_eb > 0) && (*cur_U_pos != 0) && (*cur_V_pos != 0)){
				bool unpred_flag = false;
				T decompressed[2];
				double abs_eb = log2(1 + required_eb);
				*eb_quant_index_pos = eb_exponential_quantize(abs_eb, base, log_of_base);
				// *eb_quant_index_pos = eb_linear_quantize(abs_eb, 1e-2);
				if(*eb_quant_index_pos > 0){
					// compress U and V
					for(int k=0; k<2; k++){
						T * cur_data_pos = (k == 0) ? cur_log_U_pos : cur_log_V_pos;
						T cur_data = *cur_data_pos;
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
				}
				else unpred_flag = true;
				if(unpred_flag){
					// recover quant index
					*(eb_quant_index_pos ++) = 0;
					*(data_quant_index_pos ++) = intv_radius;
					*(data_quant_index_pos ++) = intv_radius;
					unpred_data.push_back(*cur_U_pos);
					unpred_data.push_back(*cur_V_pos);
				}
				else{
					eb_quant_index_pos ++;
					data_quant_index_pos += 2;
					// assign decompressed data
					*cur_log_U_pos = decompressed[0];
					*cur_log_V_pos = decompressed[1];
					*cur_U_pos = (*cur_U_pos > 0) ? exp2(*cur_log_U_pos) : -exp2(*cur_log_U_pos);
					*cur_V_pos = (*cur_V_pos > 0) ? exp2(*cur_log_V_pos) : -exp2(*cur_log_V_pos);
				}
			}
			else{
				// record as unpredictable data
				*(eb_quant_index_pos ++) = 0;
				*(data_quant_index_pos ++) = intv_radius;
				*(data_quant_index_pos ++) = intv_radius;
				unpred_data.push_back(*cur_U_pos);
				unpred_data.push_back(*cur_V_pos);
			}
			cur_log_U_pos ++, cur_log_V_pos ++;
			cur_U_pos ++, cur_V_pos ++;
		}
	}
	// writefile("eb_2d.dat", eb, num_elements);
	// free(eb);
	free(log_U);
	free(log_V);
	free(decompressed_U);
	free(decompressed_V);
	// printf("offsets eb_q, data_q, unpred: %ld %ld %ld\n", eb_quant_index_pos - eb_quant_index, data_quant_index_pos - data_quant_index, unpred_data.size());
	unsigned char * compressed = (unsigned char *) malloc(3*num_elements*sizeof(T));
	unsigned char * compressed_pos = compressed;
	write_variable_to_dst(compressed_pos, base);
	write_variable_to_dst(compressed_pos, intv_radius);
	write_array_to_dst(compressed_pos, sign_map_compressed, 2*sign_map_size);
	free(sign_map_compressed);
	size_t unpredictable_count = unpred_data.size();
	// printf("unpredictable_count=%d\n", unpredictable_count);
	write_variable_to_dst(compressed_pos, unpredictable_count);
	write_array_to_dst(compressed_pos, (T *)&unpred_data[0], unpredictable_count);	
	Huffman_encode_tree_and_data(2*1024, eb_quant_index, num_elements, compressed_pos);
	free(eb_quant_index);
	Huffman_encode_tree_and_data(2*capacity, data_quant_index, 2*num_elements, compressed_pos);
	free(data_quant_index);
	compressed_size = compressed_pos - compressed;
	return compressed;	
}

template
unsigned char *
sz_compress_cp_preserve_2d_bilinear_online_log(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_2d_bilinear_online_log(const double * U, const double * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);