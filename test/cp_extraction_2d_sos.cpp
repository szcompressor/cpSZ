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

uint64_t vector_field_scaling_factor = 1;
int DW = 128, DH = 128;// the dimensionality of the data is DW*DH

ftk::ndarray<float> grad;
ftk::simplicial_regular_mesh m(2); // the 2D spatial mesh
std::mutex mutex;

struct critical_point_t {
  double mu[3]; // interpolation coefficients
  double x[2]; // the coordinates of the critical points
  double v;
  int type;
  size_t simplex_id;
  critical_point_t(){}
  critical_point_t(const double * mu_, const double * x_, const double v_, const int t_){
    for(int i=0; i<3; i++) mu[i] = mu_[i];
    for(int i=0; i<2; i++) x[i] = x_[i];
    v = v_;
    type = t_;
  }
};
std::unordered_map<int, critical_point_t> critical_points;
 
#define DEFAULT_EB 1
#define SINGULAR 0
#define ATTRACTING 1 // 2 real negative eigenvalues
#define REPELLING 2 // 2 real positive eigenvalues
#define SADDLE 3// 1 real negative and 1 real positive
#define ATTRACTING_FOCUS 4 // complex with negative real
#define REPELLING_FOCUS 5 // complex with positive real
#define CENTER 6 // complex with 0 real

std::string get_critical_point_type_string(int type){
  switch(type){
    case 0:
      return "SINGULAR";
    case 1:
      return "ATTRACTING";
    case 2:
      return "REPELLING";
    case 3:
      return "SADDLE";
    case 4:
      return "ATTRACTING_FOCUS";
    case 5:
      return "REPELLING_FOCUS";
    case 6:
      return "CENTER";
    default:
      return "INVALID";
  }
}

bool check_simplex(const ftk::simplicial_regular_mesh_element& s)
{
  if (!s.valid(m)) return false; // check if the 3-simplex is valid
  const auto &vertices = s.vertices(m);
  // simplex_vectors(vertices, v);
  double X[3][2], v[3][2];
  for (int i = 0; i < 3; i ++) {
    for (int j = 0; j < 2; j ++)
      v[i][j] = grad(j, vertices[i][0], vertices[i][1]);
    for (int j = 0; j < 2; j ++)
      X[i][j] = vertices[i][j];
  }

#if FTK_HAVE_GMP
  typedef mpf_class fp_t;
  fp_t vf[3][2];
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 2; j ++) {
      const double x = v[i][j];
      if (std::isnan(x) || std::isinf(x)) return false;
      else vf[i][j] = v[i][j];
    }
#else
  int64_t vf[3][2];
  for (int i = 0; i < 3; i ++)
    for (int j = 0; j < 2; j ++) {
      const double x = v[i][j];
      if (std::isnan(x) || std::isinf(x)) return false;
      else vf[i][j] = v[i][j] * vector_field_scaling_factor;
    }
#endif

  // robust critical point test
  int indices[3];
  for (int i = 0; i < vertices.size(); i ++)
    indices[i] = m.get_lattice().to_integer(vertices[i]);
  bool succ = ftk::robust_critical_point_in_simplex2(vf, indices);
  if (!succ) return false;

  double mu[3]; // check intersection
  double cond;
  bool succ2 = ftk::inverse_lerp_s2v2(v, mu, &cond);
  // if (!succ2) return false;
  // if (std::isnan(mu[0]) || std::isnan(mu[1]) || std::isnan(mu[2])) return false;
  // fprintf(stderr, "mu=%f, %f, %f\n", mu[0], mu[1], mu[2]);
  if (!succ2) ftk::clamp_barycentric<3>(mu);
  double x[2]; // position
  // simplex_coordinates(vertices, X);
  ftk::lerp_s2v2(X, mu, x);
  // fprintf(stdout, "simplex_id=%d, corner=%d, %d, type=%d, mu=%f, %f, %f, x=%f, %f\n", s.to_integer(), s.corner[0], s.corner[1], s.type, mu[0], mu[1], mu[2], x[0], x[1]);
  double J[2][2]; // jacobian
  ftk::jacobian_2dsimplex2(X, v, J);  
  int cp_type = 0;
  std::complex<double> eig[2];
  double delta = ftk::solve_eigenvalues2x2(J, eig);
  // if(fabs(delta) < std::numeric_limits<double>::epsilon())
  if (delta >= 0) { // two real roots
    if (eig[0].real() * eig[1].real() < 0) {
      cp_type = SADDLE;
    } else if (eig[0].real() < 0) {
      cp_type = ATTRACTING;
    }
    else if (eig[0].real() > 0){
      cp_type = REPELLING;
    }
    else cp_type = SINGULAR;
  } else { // two conjugate roots
    if (eig[0].real() < 0) {
      cp_type = ATTRACTING_FOCUS;
    } else if (eig[0].real() > 0) {
      cp_type = REPELLING_FOCUS;
    } else 
      cp_type = CENTER;
  }
  critical_point_t cp(mu, x, 0, cp_type);
  cp.simplex_id = s.to_integer(m);
  {
    std::lock_guard<std::mutex> guard(mutex);
    critical_points.insert(std::make_pair(s.to_integer(m), cp));
  }
  return true;
}

void extract_critical_points()
{
  fprintf(stderr, "extracting critical points...\n");
  m.set_lb_ub({1, 1}, {DW-2, DH-2}); // set the lower and upper bounds of the mesh
  m.element_for(2, check_simplex); // iterate over all 3-simplices
}

void refill_gradient(int id, const float* grad_tmp){
  const float * grad_tmp_pos = grad_tmp;
  for (int i = 0; i < DH; i ++) {
    for (int j = 0; j < DW; j ++) {
      grad(id, j, i) = *(grad_tmp_pos ++);
    }
  }
}

std::unordered_map<int, critical_point_t> get_critical_points(){
  critical_points.clear();
  // derive_gradients();
  extract_critical_points();
  std::unordered_map<int, critical_point_t> result(critical_points);
  return result;
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

void record_criticalpoints(std::string prefix, const std::vector<critical_point_t>& cps){
  double * positions = (double *) malloc(cps.size()*2*sizeof(double));
  int * type = (int *) malloc(cps.size()*sizeof(int));
  int i = 0;
  for(const auto& cp:cps){
    positions[2*i] = cp.x[0];
    positions[2*i+1] = cp.x[1];
    type[i] = cp.type;
    i ++;
  }
  writefile((prefix + "_pos.dat").c_str(), positions, cps.size()*2);
  writefile((prefix + "_type.dat").c_str(), type, cps.size());
  free(positions);
  free(type);
}

inline bool file_exists(const std::string& filename) {
    std::ifstream f(filename.c_str());
    return f.good();
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
    p.x[0] = positions[2*i]; p.x[1] = positions[2*i+1];
    p.type = type[i];
    p.simplex_id = sid[i]; 
    cps.insert(std::make_pair(sid[i], p));
  }
  return cps;
}

void record_criticalpoints(const std::string& prefix, const std::vector<critical_point_t>& cps, bool write_sid=false){
  double * positions = (double *) malloc(cps.size()*2*sizeof(double));
  int * type = (int *) malloc(cps.size()*sizeof(int));
  int i = 0;
  for(const auto& cp:cps){
    positions[2*i] = cp.x[0];
    positions[2*i+1] = cp.x[1];
    type[i] = cp.type;
    i ++;
  }
  writefile((prefix + "_pos.dat").c_str(), positions, cps.size()*2);
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

int main(int argc, char **argv){
  size_t num = 0;
  float * u = readfile<float>(argv[1], num);
  float * v = readfile<float>(argv[2], num);
  DW = atoi(argv[3]);
  DH = atoi(argv[4]);

  grad.reshape({2, static_cast<unsigned long>(DW), static_cast<unsigned long>(DH)});
  refill_gradient(0, u);
  refill_gradient(1, v);
  // compute vector_field_resolution
  const int type_bits = 63;
  double vector_field_resolution = 0;
  vector_field_scaling_factor = 1;
  for (int i=0; i<num; i++){
    double min_val = std::max(fabs(u[i]), fabs(v[i]));
    vector_field_resolution = std::max(vector_field_resolution, min_val);
  }
  int vbits = std::ceil(std::log2(vector_field_resolution));
  int nbits = (type_bits - 3) / 2;
  vector_field_scaling_factor = 1 << (nbits - vbits);
  std::cerr << "resolution=" << vector_field_resolution 
  << ", factor=" << vector_field_scaling_factor 
  << ", nbits=" << nbits << ", vbits=" << vbits << ", shift_bits=" << nbits - vbits << std::endl;

  free(u);
  free(v);

  // auto critical_points_0 = get_critical_points();
  std::string cp_prefix = "origin_M_" + std::to_string(nbits - vbits) + "_bits";
  cp_prefix += "_sos";
  bool cp_file = file_exists(cp_prefix + "_sid.dat");
  cp_file ? printf("Critical point file found!\n") : printf("Critical point Not found, recomputing\n");
  auto critical_points_0 = cp_file ? read_criticalpoints(cp_prefix) : get_critical_points();
  std::cout << "critical_points number = " << critical_points_0.size() << std::endl;
  {
    std::cout << "read decompressed data\n";
    std::string fn_u = std::string(argv[1]) + ".out";
    std::string fn_v = std::string(argv[2]) + ".out";
    u = readfile<float>(fn_u.c_str(), num);
    v = readfile<float>(fn_v.c_str(), num);
    std::cout << "read " << num << " decompressed data done\n";
  }
  refill_gradient(0, u);
  refill_gradient(1, v);
  free(u);
  free(v);
  std::cout << "refill gradients done\n";
  auto critical_points_1 = get_critical_points();
  std::cout << "decompressed critical_points number = " << critical_points_1.size() << std::endl;
  int matches = 0;
  std::vector<critical_point_t> fp, fn, ft, m;
  std::vector<critical_point_t> origin;
  // int cp_type_change = 0;
  for(const auto& p:critical_points_0){
    auto cp = p.second;
    origin.push_back(cp);
    if(critical_points_1.find(p.first) != critical_points_1.end()){
      // matches ++;
      auto cp_1 = critical_points_1[p.first];
      if(cp.type != cp_1.type){
        // std::cout << "Change from " << get_critical_point_type_string(cp.type)
        //   << " to " << get_critical_point_type_string(cp_1.type) <<
        //    ", positions from " << cp.x[0] << ", " << cp.x[1] << " to " << cp_1.x[0] << ", " << cp_1.x[1] << std::endl;
        // cp_type_change ++;
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
  // for(const auto& cp:fn){
  //   std::cout << "x = " << cp.x[0] << ", " << cp.x[1] << "; type = " << get_critical_point_type_string(cp.type) << std::endl;
  // }
  std::cout << "FP number = " << fp.size() << std::endl; 
  // for(const auto& cp:fp){
  //   std::cout << "x = " << cp.x[0] << ", " << cp.x[1] << "; type = " << get_critical_point_type_string(cp.type) << std::endl;
  // }
  std::cout << "FT number = " << ft.size() << std::endl; 
  // for(const auto& cp:ft){
  //   std::cout << "x = " << cp.x[0] << ", " << cp.x[1] << "; type = " << get_critical_point_type_string(cp.type) << std::endl;
  // }
  std::cout << "Ground truth = " << origin.size() << std::endl;
  std::cout << "TP = " << m.size() << std::endl;
  std::cout << "FP = " << fp.size() << std::endl;
  std::cout << "FN = " << fn.size() << std::endl;
  std::cout << "FT = " << ft.size() << std::endl;
  record_criticalpoints(cp_prefix, origin, true);    
  // if(!cp_file) {
  //   record_criticalpoints(std::string("origin_M"), origin, true);
  //   std::string prefix = std::string(argv[5]);
  //   record_criticalpoints(prefix+"_M", m, false);    
  //   record_criticalpoints(prefix+"_FP", fp, false);    
  //   record_criticalpoints(prefix+"_FN", fn, false);    
  //   record_criticalpoints(prefix+"_FT", ft, false);    
  // }
  return 0;
}