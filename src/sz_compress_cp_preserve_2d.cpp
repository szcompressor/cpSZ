#include "sz_cp_preserve_utils.hpp"
#include "sz_compress_3d.hpp"
#include "sz_compress_cp_preserve_2d.hpp"
#include "sz_def.hpp"
#include "sz_compression_utils.hpp"

template<typename Type>
void writefile(const char * file, Type * data, size_t num_elements){
	std::ofstream fout(file, std::ios::binary);
	fout.write(reinterpret_cast<const char*>(&data[0]), num_elements*sizeof(Type));
	fout.close();
}

// maximal error bound to keep the sign of u0v1 - u0v2 + u1v2 - u1v0 + u2v0 - u2v1
template<typename T>
inline double max_eb_to_keep_sign_det2x2(const T u0, const T u1, const T u2, const T v0, const T v1, const T v2){
  double positive = 0;
  double negative = 0;
  accumulate(u0*v1, positive, negative);
  accumulate(-u0*v2, positive, negative);
  accumulate(u1*v2, positive, negative);
  accumulate(-u1*v0, positive, negative);
  accumulate(u2*v0, positive, negative);
  accumulate(-u2*v1, positive, negative);
  return max_eb_to_keep_sign(positive, negative, 2);
}

template<typename T>
inline double max_eb_to_keep_sign_2d_offline_2(const T u0, const T u1, const int degree=2){
  double positive = 0;
  double negative = 0;
  accumulate(u0, positive, negative);
  accumulate(u1, positive, negative);
  return max_eb_to_keep_sign(positive, negative, degree);
}

template<typename T>
inline double max_eb_to_keep_sign_2d_offline_4(const T u0, const T u1, const T u2, const T u3, const int degree=2){
  double positive = 0;
  double negative = 0;
  accumulate(u0, positive, negative);
  accumulate(u1, positive, negative);
  accumulate(u2, positive, negative);
  accumulate(u3, positive, negative);
  return max_eb_to_keep_sign(positive, negative, degree);
}

// det(c) = (x0 - x2)*(y1 - y2) - (x1 - x2)*(y0 - y2)
// c0 = (y1 - y2) / det(c)   c1 = -(y0 - y2) / det(c)
// c1 = -(x1 - x2) / det(c)  c3 = (x0 - x2) / det(c)
template<typename T>
inline void get_adjugate_matrix_for_position(const T x0, const T x1, const T x2, const T y0, const T y1, const T y2, T c[4]){
  T determinant = (x0 - x2)*(y1 - y2) - (x1 - x2)*(y0 - y2);
  c[0] = (y1 - y2) / determinant;
  c[1] = -(y0 - y2) / determinant;
  c[2] = -(x1 - x2) / determinant;
  c[3] = (x0 - x2) / determinant;
  // printf("%.4g, %.2g %.2g %.2g %.2g\n", determinant, c[0], c[1], c[2], c[3]);
  // exit(0);
}

// accumulate positive and negative in (a + b + c ...)^2
template<typename T>
inline void accumulate_in_square(const std::vector<T>& coeff, double& positive, double& negative){
  for(int i=0; i<coeff.size(); i++){
    for(int j=0; j<coeff.size(); j++){
      accumulate(coeff[i]*coeff[j], positive, negative);
    }
  }
}

// maximal error bound to keep the sign of B^2 - 4C
// where  B = - (c0 * (u0 - u2) + c1 * (u1 - u2) + c2 * (v0 - v2) + c3 * (v1 - v2))
//        C = det2x2 = u0v1 - u0v2 + u1v2 - u1v0 + u2v0 - u2v1
template<typename T>
inline double max_eb_to_keep_sign_eigen_delta_2(const T u0, const T u1, const T u2, const T v0, const T v1, const T v2,
  const T x0, const T x1, const T x2, const T y0, const T y1, const T y2){
  double eb = 1;
  T c[4] = {0};
  {
    get_adjugate_matrix_for_position(x0, x1, x2, y0, y1, y2, c);
    // keep sign for B
    double positive = 0;
    double negative = 0;
    accumulate(c[0]*u0, positive, negative);
    accumulate(c[1]*u1, positive, negative);
    accumulate(-(c[0] + c[1])*u2, positive, negative);
    accumulate(c[2]*v0, positive, negative);
    accumulate(c[3]*v1, positive, negative);
    accumulate(-(c[2] + c[3])*v2, positive, negative);
    eb = max_eb_to_keep_sign(positive, negative, 1);
    // keep sign for C
    eb = MIN(eb, max_eb_to_keep_sign_det2x2(u0, u1, u2, v0, v1, v2));
  }
  T m = c[1]*c[2] - c[0]*c[3];
  T C = (-m) * (u0*v1 - u0*v2 + u1*v2 - u1*v0 + u2*v0 - u2*v1);
  if(C == 0) return 0;
  if(C < 0) return eb;
  {
    std::vector<T> coeff(6);
    coeff[0] = c[0]*u0;
    coeff[1] = c[1]*u1;
    coeff[2] = - (c[1] + c[0])*u2;
    coeff[3] = c[2]*v0;
    coeff[4] = c[3]*v1;
    coeff[5] = - (c[3] + c[2])*v2;
    // keep sign for B^2 - 4*C
    double positive = 0;
    double negative = 0;
    accumulate_in_square(coeff, positive, negative);
    accumulate(-4*m*u1*v0, positive, negative);
    accumulate(4*m*u2*v0, positive, negative);
    accumulate(4*m*u0*v1, positive, negative);
    accumulate(-4*m*u2*v1, positive, negative);
    accumulate(-4*m*u0*v2, positive, negative);
    accumulate(4*m*u1*v2, positive, negative);
    eb = MIN(eb, max_eb_to_keep_sign(positive, negative, 2));
  }
  return eb;
}

template<typename T>
inline double max_eb_to_keep_position_and_type(const T u0, const T u1, const T u2, const T v0, const T v1, const T v2,
											const T x0, const T x1, const T x2, const T y0, const T y1, const T y2){
	double u0v1 = u0 * v1;
	double u1v0 = u1 * v0;
	double u0v2 = u0 * v2;
	double u2v0 = u2 * v0;
	double u1v2 = u1 * v2;
	double u2v1 = u2 * v1;
	double det = u0v1 - u1v0 + u1v2 - u2v1 + u2v0 - u0v2;
	double eb = 0;
	if(det != 0){
		bool f1 = (det / (u2v0 - u0v2) >= 1);
		bool f2 = (det / (u1v2 - u2v1) >= 1); 
		bool f3 = (det / (u0v1 - u1v0) >= 1); 
		if(f1 && f2 && f3){
			// critical point found
			eb = 1;
			double eb1 = MIN(max_eb_to_keep_sign_2d_offline_2(u2v0, -u0v2), max_eb_to_keep_sign_2d_offline_4(u0v1, -u1v0, u1v2, -u2v1));
			double eb2 = MIN(max_eb_to_keep_sign_2d_offline_2(u1v2, -u2v1), max_eb_to_keep_sign_2d_offline_4(u0v1, -u1v0, u2v0, -u0v2));
			double eb3 = MIN(max_eb_to_keep_sign_2d_offline_2(u0v1, -u1v0), max_eb_to_keep_sign_2d_offline_4(u1v2, -u2v1, u2v0, -u0v2));
			double eb4 = MIN(eb3, max_eb_to_keep_sign_eigen_delta_2(u0, u1, u2, v0, v1, v2, x0, x1, x2, y0, y1, y2));
			eb = MIN(eb1, eb2);
			eb = MIN(eb, eb4);
		}
		else{
			// no critical point
			eb = 0;
			if(!f1){
				double eb_cur = MIN(max_eb_to_keep_sign_2d_offline_2(u2v0, -u0v2), max_eb_to_keep_sign_2d_offline_4(u0v1, -u1v0, u1v2, -u2v1));
				// double eb_cur = MIN(max_eb_to_keep_sign_2(u2, u0, v2, v0), max_eb_to_keep_sign_4(u0, u1, u2, v0, v1, v2));
				eb = MAX(eb, eb_cur);
			}
			if(!f2){
				double eb_cur = MIN(max_eb_to_keep_sign_2d_offline_2(u1v2, -u2v1), max_eb_to_keep_sign_2d_offline_4(u0v1, -u1v0, u2v0, -u0v2));
				// double eb_cur = MIN(max_eb_to_keep_sign_2(u1, u2, v1, v2), max_eb_to_keep_sign_4(u2, u0, u1, v2, v0, v1));
				eb = MAX(eb, eb_cur);
			}
			if(!f3){
				double eb_cur = MIN(max_eb_to_keep_sign_2d_offline_2(u0v1, -u1v0), max_eb_to_keep_sign_2d_offline_4(u1v2, -u2v1, u2v0, -u0v2));
				// double eb_cur = MIN(max_eb_to_keep_sign_2(u0, u1, v0, v1), max_eb_to_keep_sign_4(u1, u2, u0, v1, v2, v0));
				eb = MAX(eb, eb_cur);
			}
			// eb = MIN(eb, DEFAULT_EB);
		}
	}
	return eb;
}

// compression with pre-computed error bounds
template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_offline(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){

	size_t num_elements = r1 * r2;
	double * eb = (double *) malloc(num_elements * sizeof(double));
	for(int i=0; i<num_elements; i++) eb[i] = max_pwr_eb;
	const T * U_pos = U;
	const T * V_pos = V;
	double * eb_pos = eb;
	// coordinates for triangle_coordinates
	const T X_upper[3][2] = {{0, 0}, {1, 0}, {1, 1}};
	const T X_lower[3][2] = {{0, 0}, {0, 1}, {1, 1}};
	const size_t offset_upper[3] = {0, r2, r2+1};
	const size_t offset_lower[3] = {0, 1, r2+1};
	printf("compute eb\n");
	for(int i=0; i<r1-1; i++){
		const T * U_row_pos = U_pos;
		const T * V_row_pos = V_pos;
		double * eb_row_pos = eb_pos;
		for(int j=0; j<r2-1; j++){
			for(int k=0; k<2; k++){
				auto X = (k == 0) ? X_upper : X_lower;
				auto offset = (k == 0) ? offset_upper : offset_lower;
				// reversed order!
				double max_cur_eb = max_eb_to_keep_position_and_type(U_row_pos[offset[0]], U_row_pos[offset[1]], U_row_pos[offset[2]],
					V_row_pos[offset[0]], V_row_pos[offset[1]], V_row_pos[offset[2]], X[0][1], X[1][1], X[2][1],
					X[0][0], X[1][0], X[2][0]);
				eb_row_pos[offset[0]] = MIN(eb_row_pos[offset[0]], max_cur_eb);
				eb_row_pos[offset[1]] = MIN(eb_row_pos[offset[1]], max_cur_eb);
				eb_row_pos[offset[2]] = MIN(eb_row_pos[offset[2]], max_cur_eb);
			}
			U_row_pos ++;
			V_row_pos ++;
			eb_row_pos ++;
		}
		U_pos += r2;
		V_pos += r2;
		eb_pos += r2;
	}
	printf("compute eb done\n");
	double * eb_u = (double *) malloc(num_elements * sizeof(double));
	double * eb_v = (double *) malloc(num_elements * sizeof(double));
	int * eb_quant_index = (int *) malloc(2*num_elements*sizeof(int));
	int * eb_quant_index_pos = eb_quant_index;
	const int base = 4;
	double log2_of_base = log2(base);
	const double threshold = std::numeric_limits<double>::epsilon();
	for(int i=0; i<num_elements; i++){
		eb_u[i] = fabs(U[i]) * eb[i];
		*(eb_quant_index_pos ++) = eb_exponential_quantize(eb_u[i], base, log2_of_base, threshold);
		// *(eb_quant_index_pos ++) = eb_linear_quantize(eb_u[i], 1e-2);
		if(eb_u[i] < threshold) eb_u[i] = 0;
	}
	for(int i=0; i<num_elements; i++){
		eb_v[i] = fabs(V[i]) * eb[i];
		*(eb_quant_index_pos ++) = eb_exponential_quantize(eb_v[i], base, log2_of_base, threshold);
		// *(eb_quant_index_pos ++) = eb_linear_quantize(eb_v[i], 1e-2);
		if(eb_v[i] < threshold) eb_v[i] = 0;
	}
	free(eb);
	printf("quantize eb done\n");
	unsigned char * compressed_eb = (unsigned char *) malloc(2*num_elements*sizeof(int));
	unsigned char * compressed_eb_pos = compressed_eb; 
	Huffman_encode_tree_and_data(2*1024, eb_quant_index, 2*num_elements, compressed_eb_pos);
	size_t compressed_eb_size = compressed_eb_pos - compressed_eb;
	size_t compressed_u_size = 0;
	size_t compressed_v_size = 0;
	unsigned char * compressed_u = sz_compress_2d_with_eb(U, eb_u, r1, r2, compressed_u_size);
	unsigned char * compressed_v = sz_compress_2d_with_eb(V, eb_v, r1, r2, compressed_v_size);
	printf("eb_size = %ld, u_size = %ld, v_size = %ld\n", compressed_eb_size, compressed_u_size, compressed_v_size);
	free(eb_u);
	free(eb_v);
	compressed_size = sizeof(int) + sizeof(size_t) + compressed_eb_size + sizeof(size_t) + compressed_u_size + sizeof(size_t) + compressed_v_size;
	unsigned char * compressed = (unsigned char *) malloc(compressed_size);
	unsigned char * compressed_pos = compressed;
	write_variable_to_dst(compressed_pos, base);
	write_variable_to_dst(compressed_pos, threshold);
	write_variable_to_dst(compressed_pos, compressed_eb_size);
	write_variable_to_dst(compressed_pos, compressed_u_size);
	write_variable_to_dst(compressed_pos, compressed_v_size);
	write_array_to_dst(compressed_pos, compressed_eb, compressed_eb_size);
	write_array_to_dst(compressed_pos, compressed_u, compressed_u_size);
	printf("compressed_pos = %ld\n", compressed_pos - compressed);
	write_array_to_dst(compressed_pos, compressed_v, compressed_v_size);
	free(compressed_eb);
	free(compressed_u);
	free(compressed_v);
	return compressed;
}

template
unsigned char *
sz_compress_cp_preserve_2d_offline(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_2d_offline(const double * U, const double * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

// compression with pre-computed error bounds in logarithmic domain
template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_offline_log(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){

	size_t num_elements = r1 * r2;
	double * eb = (double *) malloc(num_elements * sizeof(double));
	for(int i=0; i<num_elements; i++) eb[i] = max_pwr_eb;
	const T * U_pos = U;
	const T * V_pos = V;
	double * eb_pos = eb;
	// coordinates for triangle_coordinates
	const T X_upper[3][2] = {{0, 0}, {1, 0}, {1, 1}};
	const T X_lower[3][2] = {{0, 0}, {0, 1}, {1, 1}};
	const size_t offset_upper[3] = {0, r2, r2+1};
	const size_t offset_lower[3] = {0, 1, r2+1};
	// printf("compute eb\n");
	for(int i=0; i<r1-1; i++){
		const T * U_row_pos = U_pos;
		const T * V_row_pos = V_pos;
		double * eb_row_pos = eb_pos;
		for(int j=0; j<r2-1; j++){
			for(int k=0; k<2; k++){
				auto X = (k == 0) ? X_upper : X_lower;
				auto offset = (k == 0) ? offset_upper : offset_lower;
				// reversed order!
				double max_cur_eb = max_eb_to_keep_position_and_type(U_row_pos[offset[0]], U_row_pos[offset[1]], U_row_pos[offset[2]],
					V_row_pos[offset[0]], V_row_pos[offset[1]], V_row_pos[offset[2]], X[0][1], X[1][1], X[2][1],
					X[0][0], X[1][0], X[2][0]);
				eb_row_pos[offset[0]] = MIN(eb_row_pos[offset[0]], max_cur_eb);
				eb_row_pos[offset[1]] = MIN(eb_row_pos[offset[1]], max_cur_eb);
				eb_row_pos[offset[2]] = MIN(eb_row_pos[offset[2]], max_cur_eb);
			}
			U_row_pos ++;
			V_row_pos ++;
			eb_row_pos ++;
		}
		U_pos += r2;
		V_pos += r2;
		eb_pos += r2;
	}
	// writefile("eb_2d_decoupled.dat", eb, num_elements);
	// printf("compute eb done\n");
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
	// transfrom eb to log(1 + eb) and the quantize
	int * eb_quant_index = (int *) malloc(num_elements*sizeof(int));
	int * eb_quant_index_pos = eb_quant_index;
	const int base = 2;
	double log2_of_base = log2(base);
	const double threshold = std::numeric_limits<double>::epsilon();
	for(int i=0; i<num_elements; i++){
		eb[i] = log2(1 + eb[i]);
		*(eb_quant_index_pos ++) = eb_exponential_quantize(eb[i], base, log2_of_base, threshold);
		// *(eb_quant_index_pos ++) = eb_linear_quantize(eb[i], 5e-3);
	}
	// printf("quantize eb done\n");
	unsigned char * compressed_eb = (unsigned char *) malloc(num_elements*sizeof(int));
	unsigned char * compressed_eb_pos = compressed_eb; 
	Huffman_encode_tree_and_data(2*1024, eb_quant_index, num_elements, compressed_eb_pos);
	size_t compressed_eb_size = compressed_eb_pos - compressed_eb;
	size_t compressed_u_size = 0;
	size_t compressed_v_size = 0;
	unsigned char * compressed_u = sz_compress_2d_with_eb(log_U, eb, r1, r2, compressed_u_size);
	free(log_U);
	unsigned char * compressed_v = sz_compress_2d_with_eb(log_V, eb, r1, r2, compressed_v_size);
	free(log_V);
	// printf("eb_size = %ld, log_u_size = %ld, log_v_size = %ld\n", compressed_eb_size, compressed_u_size, compressed_v_size);
	free(eb);
	compressed_size = sizeof(int) + 2*sign_map_size + sizeof(size_t) + sizeof(double) + compressed_eb_size + sizeof(size_t) + compressed_u_size + sizeof(size_t) + compressed_v_size;
	unsigned char * compressed = (unsigned char *) malloc(compressed_size);
	unsigned char * compressed_pos = compressed;
	write_variable_to_dst(compressed_pos, base);
	write_variable_to_dst(compressed_pos, threshold);
	write_variable_to_dst(compressed_pos, compressed_eb_size);
	write_variable_to_dst(compressed_pos, compressed_u_size);
	write_variable_to_dst(compressed_pos, compressed_v_size);
	write_array_to_dst(compressed_pos, compressed_eb, compressed_eb_size);
	write_array_to_dst(compressed_pos, sign_map_compressed, 2*sign_map_size);
	// printf("before data: %ld\n", compressed_pos - compressed);
	write_array_to_dst(compressed_pos, compressed_u, compressed_u_size);
	write_array_to_dst(compressed_pos, compressed_v, compressed_v_size);
	free(sign_map_compressed);
	free(compressed_eb);
	free(compressed_u);
	free(compressed_v);
	return compressed;
}

template
unsigned char *
sz_compress_cp_preserve_2d_offline_log(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_2d_offline_log(const double * U, const double * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

typedef struct conditions_2d{
	bool computed;
	bool singular;
	bool flags[3];
	conditions_2d(){
		computed = false;
		singular = false;
		for(int i=0; i<3; i++){
			flags[i] = false;
		}
	}
}conditions_2d;

// maximal error bound to keep the sign of A*(1 + e_1) + B*(1 + e_2) + C
template<typename T>
inline double max_eb_to_keep_sign_2d_online(const T A, const T B, const T C=0){
	double fabs_sum = (fabs(A) + fabs(B));
	if(fabs_sum == 0) return 0;
	return fabs(A + B + C) / fabs_sum;
}

// W0 + W1 = u1v2 - u2v1 + u2v0 - u0v2
// W1 + W2 = u2v0 - u0v2 + u0v1 - u1v0
// W2 + W0 = u0v1 - u1v0 + u1v2 - u2v1
template<typename T>
inline double max_eb_to_keep_position_online(const T u0v1, const T u1v0, const T u1v2, const T u2v1, const T u2v0, const T u0v2){
	double eb = MIN(max_eb_to_keep_sign_2d_online(-u2v1, u1v2), max_eb_to_keep_sign_2d_online(u2v0, -u0v2));
	// eb = MIN(eb, max_eb_to_keep_sign_2d_online(u2v0, -u0v2, u0v1 - u1v0));
	// eb = MIN(eb, max_eb_to_keep_sign_2d_online(-u2v1, u1v2, u0v1 - u1v0));
	// eb = MIN(eb, max_eb_to_keep_sign_2d_online(u2v0 - u2v1, u1v2 - u0v2));
	return eb;
}

// maximal error bound to keep the sign of B^2 - 4C
// where  B = - (c0 * (u0 - u2) + c1 * (u1 - u2) + c2 * (v0 - v2) + c3 * (v1 - v2))
//        C = det2x2 = u0v1 - u0v2 + u1v2 - u1v0 + u2v0 - u2v1
template<typename T>
inline double max_eb_to_keep_type_online(const T u0, const T u1, const T u2, const T v0, const T v1, const T v2, const T c[4]){
	double eb = 1;
	{
		// keep sign for B
	    // coeff[0] = c[0]*u0;
	    // coeff[1] = c[1]*u1;
	    // coeff[2] = - (c[1] + c[0])*u2;
	    // coeff[3] = c[2]*v0;
	    // coeff[4] = c[3]*v1;
	    // coeff[5] = - (c[3] + c[2])*v2;
		eb = max_eb_to_keep_sign_2d_online(-c[0]*u2 - c[1]*u2, -c[2]*v2 - c[3]*v2, c[0]*u0 + c[1]*u1 + c[2]*v0 + c[3]*v1);
		// keep sign for C
		eb = MIN(eb, max_eb_to_keep_sign_2d_online(u2*v0 - u2*v1, u1*v2 - u0*v2, u0*v1 - u1*v0));
	}
	T m = c[1]*c[2] - c[0]*c[3];
	T C = (-m) * (u0*v1 - u0*v2 + u1*v2 - u1*v0 + u2*v0 - u2*v1);
	if(C <= 0) return eb;
	{
		// Note that meaning of B in the rhs changes here
		// keep sign for B^2 - 4*C
		// B = A*(1+e_1) + B*(1+e_2) + C
		// C = D*(1+e_1) + E*(1+e_2) + F
		double A = -c[0]*u2 - c[1]*u2, B = -c[2]*v2 - c[3]*v2, C = c[0]*u0 + c[1]*u1 + c[2]*v0 + c[3]*v1;
		double D = (-m)*(u2*v0 - u2*v1), E = (-m)*(u1*v2 - u0*v2), F = (-m)*(u0*v1 - u1*v0);
		// B = A*e_1 + B*e_2 + C'
		// C = D*e_1 + E*e_2 + F'
		C += A + B, F += D + E;
		// B^2 - 4C = (A*e_1 + B*e_2)^2 + (2AC' - 4D)e_1 + (2BC' - 4E)e_2 + C'^2 - 4F'
		double delta = C*C - 4*F;
		if(delta == 0) return 0;
		else if(delta > 0){
			// (|2AC' - 4D| + |2BC' - 4E|)* -e + delta > 0
			if((fabs(2*A*C - 4*D) + fabs(2*B*C - 4*E)) == 0) eb = 1;
			else eb = MIN(eb, delta/(fabs(2*A*C - 4*D) + fabs(2*B*C - 4*E)));
			// check four edges

		}
		else{
			// (|A| + |B|)^2*e^2 + (|2AC' - 4D| + |2BC' - 4E|)*e + delta < 0
			// double a = (fabs(A) + fabs(B))*(fabs(A) + fabs(B));
			// double b = fabs(2*A*C - 4*D) + fabs(2*B*C - 4*E);
			// double c = delta;
			// if(b*b - 4*a*c < 0){
			// 	printf("impossible as a*c is always less than 0\n");
			// 	exit(0);
			// }
			// eb = MIN(eb, (-b + sqrt(b*b - 4*a*c))/(2*a));

			// check four corners
			double e1[2] = {-1, 1};
			double e2[2] = {-1, 1};
			double c = delta;
			for(int i=0; i<2; i++){
				for(int j=0; j<2; j++){
					double a = (e1[i] * A + e2[j] * B) * (e1[i] * A + e2[j] * B);
					double b = (2*A*C - 4*D) * e1[i] + (2*B*C - 4*E) * e2[j];
					if(a == 0) eb = MIN(eb, 1);
					else eb = MIN(eb, (-b + sqrt(b*b - 4*a*c))/(2*a));
				}
			}
		}
	}
	return eb;
}

/*
triangle mesh x0, x1, x2, derive cp-preserving eb for x2 given x0, x1
*/
template<typename T>
double 
derive_cp_eb_for_positions_online(const T u0, const T u1, const T u2, const T v0, const T v1, const T v2, const T c[4]){//, conditions_2d& cond){
	// if(!cond.computed){
	//     double M0 = u2*v0 - u0*v2;
	//     double M1 = u1*v2 - u2*v1;
	//     double M2 = u0*v1 - u1*v0;
	//     double M = M0 + M1 + M2;
	//     cond.singular = (M == 0);
	//     if(cond.singular) return 0;
	//     cond.flags[0] = (M0 == 0) || (M / M0 >= 1);
	//     cond.flags[1] = (M1 == 0) || (M / M1 >= 1);
	//     cond.flags[2] = (M2 == 0) || (M / M2 >= 1);
	//     cond.computed = true;
	// }
	// else{
	//     if(cond.singular) return 0;
	// }
	// const bool * flag = cond.flags;
	// bool f1 = flag[0];
	// bool f2 = flag[1]; 
	// bool f3 = flag[2];
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
		// eb = max_eb_to_keep_position_online(u0v1, u1v0, u1v2, u2v1, u2v0, u0v2);
		eb = MIN(max_eb_to_keep_position_online(u0*v1, u1*v0, u1*v2, u2*v1, u2*v0, u0*v2), 
			max_eb_to_keep_type_online(u0, u1, u2, v0, v1, v2, c));
	}
	else{
		eb = 0;
		if(!f1){
			// W1(W0 + W2)
			double cur_eb = MIN(max_eb_to_keep_sign_2d_online(u2*v0, -u0*v2), max_eb_to_keep_sign_2d_online(-u2*v1, u1*v2, u0*v1 - u1*v0));
			eb = MAX(eb, cur_eb);
		}
		if(!f2){
			// W0(W1 + W2)
			double cur_eb = MIN(max_eb_to_keep_sign_2d_online(-u2*v1, u1*v2), max_eb_to_keep_sign_2d_online(u2*v0, -u0*v2, u0*v1 - u1*v0));
			eb = MAX(eb, cur_eb);				
		}
		if(!f3){
			// W2(W0 + W1)
			double cur_eb = max_eb_to_keep_sign_2d_online(u2*v0 - u2*v1, u1*v2 - u0*v2);
			eb = MAX(eb, cur_eb);				
		}
	}
	return eb;
}

template<typename T>
inline bool 
inbound(T index, T lb, T ub){
	return (index >= lb) && (index < ub);
}

template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_online(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){
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
	const int base = 4;
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
	T inv_C[6][4];
	for(int i=0; i<6; i++){
		get_adjugate_matrix_for_position(x[i][0], x[i][1], x[i][2], y[i][0], y[i][1], y[i][2], inv_C[i]);
	}
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
					required_eb = MIN(required_eb, derive_cp_eb_for_positions_online(cur_U_pos[offsets[k]], cur_U_pos[offsets[k+1]], cur_U_pos[0],
						cur_V_pos[offsets[k]], cur_V_pos[offsets[k+1]], cur_V_pos[0], inv_C[k]));
				}
			}
			if(required_eb > 0){
				bool unpred_flag = false;
				T decompressed[2];
				// compress U and V
				for(int k=0; k<2; k++){
					T * cur_data_pos = (k == 0) ? cur_U_pos : cur_V_pos;
					T cur_data = *cur_data_pos;
					double abs_eb = fabs(cur_data) * required_eb;
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
sz_compress_cp_preserve_2d_online(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_2d_online(const double * U, const double * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_online_log(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){

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
	// offsets to get six adjacent triangle indices
	// the 7-th rolls back to T0
	/*
	|		T3	T4
	y	T2	X 	T5
	|	T1	T0(T6)
		-	x 	-
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
	T inv_C[6][4];
	for(int i=0; i<6; i++){
		get_adjugate_matrix_for_position(x[i][0], x[i][1], x[i][2], y[i][0], y[i][1], y[i][2], inv_C[i]);
	}
	int index_offset[6][2][2];
	for(int i=0; i<6; i++){
		for(int j=0; j<2; j++){
			index_offset[i][j][0] = x[i][j] - x[i][2];
			index_offset[i][j][1] = y[i][j] - y[i][2];
		}
	}
	// int simplex_offset[6];
	// {
	// 	// lower simplex -> 0, upper simplex -> 1
	// 	simplex_offset[0] = - 2*r2 - 2; 
	// 	simplex_offset[1] = - 2*r2 - 2 + 1; 
	// 	simplex_offset[2] = - 2; 
	// 	simplex_offset[3] = 1; 
	// 	simplex_offset[4] = 0; 
	// 	simplex_offset[5] = -2*r2 + 1; 
	// }
	// conditions_2d * conds = (conditions_2d *) malloc(2*num_elements * sizeof(conditions_2d));
	// for(int i=0; i<2*num_elements; i++) conds[i].computed = false;
	// int * count = (int *) malloc(2*num_elements*sizeof(int));
	// memset(count, 0, 2*num_elements*sizeof(int));
	// double * eb = (double *) malloc(num_elements*sizeof(double));
	// for(int i=0; i<num_elements; i++) eb[i] = 1;
	// int index_eb = 0;
	T * cur_log_U_pos = log_U;
	T * cur_log_V_pos = log_V;
	T * cur_U_pos = decompressed_U;
	T * cur_V_pos = decompressed_V;
	for(int i=0; i<r1; i++){
		// printf("start %d row\n", i);
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
					// int index = simplex_offset[k] + 2*i*r2 + 2*j;
					// conds[index].computed = false;
					// count[index] ++;
					required_eb = MIN(required_eb, derive_cp_eb_for_positions_online(cur_U_pos[offsets[k]], cur_U_pos[offsets[k+1]], cur_U_pos[0],
						cur_V_pos[offsets[k]], cur_V_pos[offsets[k+1]], cur_V_pos[0], inv_C[k]));
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
	// printf("%d %d\n", index_eb, num_elements);
	// writefile("eb_2d.dat", eb, num_elements);
	// free(eb);
	free(log_U);
	free(log_V);
	free(decompressed_U);
	free(decompressed_V);
	// printf("offsets eb_q, data_q, unpred: %ld %ld %ld\n", eb_quant_index_pos - eb_quant_index, data_quant_index_pos - data_quant_index, unpred_data.size());
	unsigned char * compressed = (unsigned char *) malloc(2*num_elements*sizeof(T));
	unsigned char * compressed_pos = compressed;
	write_variable_to_dst(compressed_pos, base);
	write_variable_to_dst(compressed_pos, intv_radius);
	write_array_to_dst(compressed_pos, sign_map_compressed, 2*sign_map_size);
	free(sign_map_compressed);
	size_t unpredictable_count = unpred_data.size();
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
sz_compress_cp_preserve_2d_online_log(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_2d_online_log(const double * U, const double * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);
