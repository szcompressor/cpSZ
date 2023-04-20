#include <mutex>
#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/bilinear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/algorithms/cca.hh>
#include <ftk/geometry/cc2curves.hh>
#include <ftk/geometry/curve2tube.hh>
#include <ftk/ndarray.hh>
#include <ftk/mesh/simplicial_regular_mesh.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <unordered_map>
#include <vector>
#include <queue>
#include <fstream>
#include <ctime>

struct critical_point_t {
  double x[3];
  int type;
  size_t simplex_id;
  critical_point_t(){}
};
 
#define SINGULAR 0
#define STABLE_SOURCE 1
#define UNSTABLE_SOURCE 2
#define STABLE_REPELLING_SADDLE 3
#define UNSTABLE_REPELLING_SADDLE 4
#define STABLE_ATRACTTING_SADDLE  5
#define UNSTABLE_ATRACTTING_SADDLE  6
#define STABLE_SINK 7
#define UNSTABLE_SINK 8

static const int tet_coords[6][4][3] = {
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

template<typename T_acc, typename T>
static inline void 
update_value(T_acc v[4][3], int local_id, int global_id, const T * U, const T * V, const T * W){
  v[local_id][0] = U[global_id];
  v[local_id][1] = V[global_id];
  v[local_id][2] = W[global_id];
}

template<typename T, typename T_fp>
static inline void 
update_index_and_value(double v[4][3], int64_t vf[4][3], int indices[4], int local_id, int global_id, const T * U, const T * V, const T * W, const T_fp * U_fp, const T_fp * V_fp, const T_fp * W_fp){
  indices[local_id] = global_id;
  update_value(v, local_id, global_id, U, V, W);
  update_value(vf, local_id, global_id, U_fp, V_fp, W_fp);
}

template<typename T_fp>
static void 
check_simplex_seq(const T_fp vf[4][3], const double v[4][3], const double X[3][3], const int indices[4], int i, int j, int k, int simplex_id, std::unordered_map<int, critical_point_t>& critical_points){
  // robust critical point test
  bool succ = ftk::robust_critical_point_in_simplex3(vf, indices);
  if (!succ) return;
  double mu[4]; // check intersection
  double cond;
  bool succ2 = ftk::inverse_lerp_s3v3(v, mu, &cond);
  if(!succ2) ftk::clamp_barycentric<4>(mu);
  double x[3]; // position
  ftk::lerp_s3v3(X, mu, x);
  critical_point_t cp;
  cp.x[0] = k + x[0]; cp.x[1] = j + x[1]; cp.x[2] = i + x[2];
  cp.type = get_cp_type(X, v);
  cp.simplex_id = simplex_id;
  critical_points[simplex_id] = cp;
}

template<typename T>
std::unordered_map<int, critical_point_t>
compute_critical_points(const T * U, const T * V, const T * W, int r1, int r2, int r3, uint64_t vector_field_scaling_factor){
  // check cp for all cells
  using T_fp = int64_t;
  ptrdiff_t dim0_offset = r2*r3;
  ptrdiff_t dim1_offset = r3;
  ptrdiff_t cell_dim0_offset = (r2-1)*(r3-1);
  ptrdiff_t cell_dim1_offset = r3-1;
  size_t num_elements = r1*r2*r3;
  T_fp * U_fp = (T_fp *) malloc(num_elements * sizeof(T_fp));
  T_fp * V_fp = (T_fp *) malloc(num_elements * sizeof(T_fp));
  T_fp * W_fp = (T_fp *) malloc(num_elements * sizeof(T_fp));
  for(int i=0; i<num_elements; i++){
    U_fp[i] = U[i]*vector_field_scaling_factor;
    V_fp[i] = V[i]*vector_field_scaling_factor;
    W_fp[i] = W[i]*vector_field_scaling_factor;
  }
  int indices[4] = {0};
  // __int128 vf[4][3] = {0};
  int64_t vf[4][3] = {0};
  double v[4][3] = {0};
  double actual_coords[6][4][3];
  for(int i=0; i<6; i++){
    for(int j=0; j<4; j++){
      for(int k=0; k<3; k++){
        actual_coords[i][j][k] = tet_coords[i][j][k];
      }
    }
  }
  std::unordered_map<int, critical_point_t> critical_points;
  for(int i=0; i<r1-1; i++){
    if(i%10==0) std::cout << i << " / " << r1-1 << std::endl;
    for(int j=0; j<r2-1; j++){
      for(int k=0; k<r3-1; k++){
        // order (reserved, z->x):
        // ptrdiff_t cell_offset = 6*(i*cell_dim0_offset + j*cell_dim1_offset + k);
        // ftk index
        ptrdiff_t cell_offset = 6*(i*dim0_offset + j*dim1_offset + k);
        // (ftk-0) 000, 001, 011, 111
        update_index_and_value(v, vf, indices, 0, i*dim0_offset + j*dim1_offset + k, U, V, W, U_fp, V_fp, W_fp);
        update_index_and_value(v, vf, indices, 1, (i+1)*dim0_offset + j*dim1_offset + k, U, V, W, U_fp, V_fp, W_fp);
        update_index_and_value(v, vf, indices, 2, (i+1)*dim0_offset + (j+1)*dim1_offset + k, U, V, W, U_fp, V_fp, W_fp);
        update_index_and_value(v, vf, indices, 3, (i+1)*dim0_offset + (j+1)*dim1_offset + (k+1), U, V, W, U_fp, V_fp, W_fp);
        check_simplex_seq(vf, v, actual_coords[0], indices, i, j, k, cell_offset, critical_points);
        // (ftk-2) 000, 010, 011, 111
        update_index_and_value(v, vf, indices, 1, i*dim0_offset + (j+1)*dim1_offset + k, U, V, W, U_fp, V_fp, W_fp);
        check_simplex_seq(vf, v, actual_coords[1], indices, i, j, k, cell_offset + 2, critical_points);
        // (ftk-1) 000, 001, 101, 111
        update_index_and_value(v, vf, indices, 1, (i+1)*dim0_offset + j*dim1_offset + k, U, V, W, U_fp, V_fp, W_fp);
        update_index_and_value(v, vf, indices, 2, (i+1)*dim0_offset + j*dim1_offset + k+1, U, V, W, U_fp, V_fp, W_fp);
        check_simplex_seq(vf, v, actual_coords[2], indices, i, j, k, cell_offset + 1, critical_points);
        // (ftk-4) 000, 100, 101, 111
        update_index_and_value(v, vf, indices, 1, i*dim0_offset + j*dim1_offset + k+1, U, V, W, U_fp, V_fp, W_fp);
        check_simplex_seq(vf, v, actual_coords[3], indices, i, j, k, cell_offset + 4, critical_points);
        // (ftk-3) 000, 010, 110, 111
        update_index_and_value(v, vf, indices, 1, i*dim0_offset + (j+1)*dim1_offset + k, U, V, W, U_fp, V_fp, W_fp);
        update_index_and_value(v, vf, indices, 2, i*dim0_offset + (j+1)*dim1_offset + k+1, U, V, W, U_fp, V_fp, W_fp);
        check_simplex_seq(vf, v, actual_coords[4], indices, i, j, k, cell_offset + 3, critical_points);
        // (ftk-5) 000, 100, 110, 111
        update_index_and_value(v, vf, indices, 1, i*dim0_offset + j*dim1_offset + k+1, U, V, W, U_fp, V_fp, W_fp);
        check_simplex_seq(vf, v, actual_coords[5], indices, i, j, k, cell_offset + 5, critical_points);
      }
    }
  }
  free(U_fp);
  free(V_fp);
  free(W_fp);
  return critical_points; 
}

template<typename Type>
Type * readfile(const char * file, size_t& num){
  std::ifstream fin(file, std::ios::binary);
  if(!fin){
        std::cout << " Error, Couldn't find the file" << "\n";
        return 0;
    }
    fin.seekg(0, std::ios::end);
    const size_t num_elements = fin.tellg() / sizeof(Type);
    fin.seekg(0, std::ios::beg);
    Type * data = (Type *) malloc(num_elements*sizeof(Type));
  fin.read(reinterpret_cast<char*>(&data[0]), num_elements*sizeof(Type));
  fin.close();
  num = num_elements;
  return data;
}

template<typename Type>
void writefile(const char * file, Type * data, size_t num_elements){
  std::ofstream fout(file, std::ios::binary);
  fout.write(reinterpret_cast<const char*>(&data[0]), num_elements*sizeof(Type));
  fout.close();
}

void record_criticalpoints(const std::string& prefix, const std::vector<critical_point_t>& cps, bool write_sid=false){
  double * positions = (double *) malloc(cps.size()*3*sizeof(double));
  int * type = (int *) malloc(cps.size()*sizeof(int));
  int i = 0;
  for(const auto& cp:cps){
    positions[3*i] = cp.x[0];
    positions[3*i+1] = cp.x[1];
    positions[3*i+2] = cp.x[2];
    type[i] = cp.type;
    i ++;
  }
  writefile((prefix + "_pos.dat").c_str(), positions, cps.size()*3);
  writefile((prefix + "_type.dat").c_str(), type, cps.size());
  if(write_sid){
    size_t * sid = (size_t *) malloc(cps.size()*sizeof(size_t));
    int i = 0;
    for(const auto& cp:cps){
      sid[i ++] = cp.simplex_id;
    }
    writefile((prefix + "_sid.dat").c_str(), sid, cps.size());
    free(sid);
  }
  free(positions);
  free(type);
}

std::unordered_map<int, critical_point_t> read_criticalpoints(const std::string& prefix){
  std::unordered_map<int, critical_point_t> cps;
  size_t num = 0;
  double * positions = readfile<double>((prefix + "_pos.dat").c_str(), num);
  int * type = readfile<int>((prefix + "_type.dat").c_str(), num);
  size_t * sid = readfile<size_t>((prefix + "_sid.dat").c_str(), num);
  printf("Read %ld critical points\n", num);
  for(int i=0; i<num; i++){
    critical_point_t p;
    p.x[0] = positions[3*i]; p.x[1] = positions[3*i+1]; p.x[2] = positions[3*i+2];
    p.type = type[i];
    p.simplex_id = sid[i]; 
    cps.insert(std::make_pair(sid[i], p));
  }
  return cps;
}

inline bool file_exists(const std::string& filename) {
    std::ifstream f(filename.c_str());
    return f.good();
}

int main(int argc, char **argv)
{
  size_t num = 0;
  fprintf(stderr, "reading inputs...\n");
  float * u = readfile<float>(argv[1], num);
  float * v = readfile<float>(argv[2], num);
  float * w = readfile<float>(argv[3], num);

  int DW = atoi(argv[4]);
  int DH = atoi(argv[5]);
  int DD = atoi(argv[6]);
  int sos = 1;
  // compute vector_field_resolution
  const int type_bits = 63;
  double vector_field_resolution = 0;
  uint64_t vector_field_scaling_factor = 1;
  for (int i=0; i<num; i++){
    double min_val = std::max(std::max(fabs(u[i]), fabs(v[i])), fabs(w[i]));
    vector_field_resolution = std::max(vector_field_resolution, min_val);
  }
  int vbits = std::ceil(std::log2(vector_field_resolution));
  int nbits = (type_bits - 5) / 3;
  vector_field_scaling_factor = 1 << (nbits - vbits);
  std::cerr << "resolution=" << vector_field_resolution 
  << ", factor=" << vector_field_scaling_factor 
  << ", nbits=" << nbits << ", vbits=" << vbits << ", shift_bits=" << nbits - vbits << std::endl;

  std::string cp_prefix = "origin_" + std::to_string(nbits - vbits) + "_bits";
  if(sos) cp_prefix += "_sos_M";
  bool cp_file = file_exists(cp_prefix + "_sid.dat");
  // cp_file = false;
  cp_file ? printf("Critical point file found!\n") : printf("Critical point Not found, recomputing\n");
  auto critical_points_0 = cp_file ? read_criticalpoints(cp_prefix) : compute_critical_points(u, v, w, DD, DH, DW, vector_field_scaling_factor);

  free(u);
  free(v);
  free(w);

  fprintf(stderr, "reading decompressed data...\n");
  std::string fn_u = std::string(argv[1]) + ".out";
  std::string fn_v = std::string(argv[2]) + ".out";
  std::string fn_w = std::string(argv[3]) + ".out";

  u = readfile<float>(fn_u.c_str(), num);
  v = readfile<float>(fn_v.c_str(), num);
  w = readfile<float>(fn_w.c_str(), num);

  struct timespec start, end;
  int err = 0;
  err = clock_gettime(CLOCK_REALTIME, &start);
  auto critical_points_1 = compute_critical_points(u, v, w, DD, DH, DW, vector_field_scaling_factor);
  err = clock_gettime(CLOCK_REALTIME, &end);
  std::cout << "extraction time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << std::endl;
  fprintf(stderr, "#critical_points = %lu\n", critical_points_0.size());
  fprintf(stderr, "#decompressed_critical_points = %lu\n", critical_points_1.size());

  int matches = 0;
  std::vector<critical_point_t> fp, fn, ft, m;
  std::vector<critical_point_t> origin;
  for(const auto& p:critical_points_0){
    auto cp = p.second;
    origin.push_back(p.second);
    if(critical_points_1.find(p.first) != critical_points_1.end()){
      matches ++;
      auto cp_1 = critical_points_1[p.first];
      // std::cout << "critical points in cell " << p.first << ": positions from (" << cp.x[0] << ", " << cp.x[1] << ", " << cp.x[2] << ") to (" << cp_1.x[0] << ", " << cp_1.x[1] << ", " << cp_1.x[2] << ")" << std::endl;
      if(cp.type != cp_1.type){
        // std::cout << "Change from " << cp.type
        //   << " to " << cp_1.type <<
        //    ", positions from " << cp.x[0] << ", " << cp.x[1] << ", " << cp.x[2] << " to " << cp_1.x[0] << ", " << cp_1.x[1] << ", " << cp_1.x[2] << std::endl;
        // printf("Type change, position from %.6f %.6f %.6f to %.6f %.6f %.6f\n", cp.x[0], cp.x[1], cp.x[2], cp_1.x[0], cp_1.x[1], cp_1.x[2]);
        ft.push_back(cp_1);
      }
      else m.push_back(cp_1);
    }
    else fn.push_back(p.second);
    // std::cout << std::endl;
  }
  for(const auto& p:critical_points_1){
    if(critical_points_0.find(p.first) == critical_points_0.end()){
      fp.push_back(p.second);
    }
  }
  std::cout << "FP number = " << fp.size() << std::endl; 
  for(const auto& cp:fp){
    // std::cout << "mu = " << cp.mu[0] << ", " << cp.mu[1] << ", " << cp.mu[2] << "; x = " << cp.x[0] << ", " << cp.x[1] << "; v = " << cp.v << std::endl;
    // std::cout << "x = " << cp.x[0] << ", " << cp.x[1] << " " << cp.x[2] << std::endl;
    printf("simplex_id = %zu, x = %.6f %.6f %.6f\n", cp.simplex_id, cp.x[0], cp.x[1], cp.x[2]);
  }
  std::cout << "FN number = " << fn.size() << std::endl; 
  for(const auto& cp:fn){
    // std::cout << "mu = " << cp.mu[0] << ", " << cp.mu[1] << ", " << cp.mu[2] << "; x = " << cp.x[0] << ", " << cp.x[1] << "; v = " << cp.v << std::endl;
    // std::cout << "x = " << cp.x[0] << ", " << cp.x[1] << " " << cp.x[2] << std::endl;
    printf("simplex_id = %zu, x = %.6f %.6f %.6f\n", cp.simplex_id, cp.x[0], cp.x[1], cp.x[2]);
  }
  std::cout << "FT number = " << ft.size() << std::endl; 
  for(const auto& cp:ft){
    // std::cout << "mu = " << cp.mu[0] << ", " << cp.mu[1] << ", " << cp.mu[2] << "; x = " << cp.x[0] << ", " << cp.x[1] << "; v = " << cp.v << std::endl;
    // std::cout << "x = " << cp.x[0] << ", " << cp.x[1] << " " << cp.x[2] << std::endl;
    printf("simplex_id = %zu, x = %.6f %.6f %.6f\n", cp.simplex_id, cp.x[0], cp.x[1], cp.x[2]);
  }
  std::cout << "Ground truth = " << critical_points_0.size() << std::endl;
  std::cout << "TP = " << m.size() << std::endl;
  std::cout << "FP = " << fp.size() << std::endl;
  std::cout << "FN = " << fn.size() << std::endl;
  std::cout << "FT = " << ft.size() << std::endl;
  // {
  // if(!cp_file) record_criticalpoints(cp_prefix, origin, true);    
  // std::string prefix = std::string(argv[7]);
  // if(m.size()) record_criticalpoints(prefix+"_M", m);    
  // if(fp.size()) record_criticalpoints(prefix+"_FP", fp);    
  // if(fn.size()) record_criticalpoints(prefix+"_FN", fn);    
  // if(ft.size()) record_criticalpoints(prefix+"_FT", ft);    
  // }
  free(u);
  free(v);
  free(w);
  return 0;
}
