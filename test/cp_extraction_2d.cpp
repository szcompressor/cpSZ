#include <ftk/numeric/print.hh>
#include <ftk/numeric/cross_product.hh>
#include <ftk/numeric/vector_norm.hh>
#include <ftk/numeric/linear_interpolation.hh>
#include <ftk/numeric/inverse_linear_interpolation_solver.hh>
#include <ftk/numeric/inverse_bilinear_interpolation_solver.hh>
#include <ftk/numeric/gradient.hh>
#include <ftk/numeric/critical_point_type.hh>
#include <ftk/numeric/critical_point_test.hh>
#include <ftk/numeric/clamp.hh>
#include <unordered_map>
#include <vector>
#include <queue>
#include <fstream>

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

static void 
  check_simplex_seq(const double v[3][2], const double X[3][2], const int indices[3], int i, int j, int simplex_id, std::unordered_map<int, critical_point_t>& critical_points){
  int zero_count = 0;
  for(int i=0; i<3; i++){
    if((v[i][0] == 0) && (v[i][1] == 0)){
      zero_count ++;
    }
  }
  if(zero_count > 0) return;
  double mu[3]; // check intersection
  double cond;
  bool succ2 = ftk::inverse_lerp_s2v2(v, mu, &cond, 0.0);
  if (!succ2) return;
  double x[2]; // position
  ftk::lerp_s2v2(X, mu, x);
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
  critical_point_t cp;
  cp.x[0] = j + x[0]; cp.x[1] = i + x[1];
  cp.simplex_id = simplex_id;
  cp.type = cp_type;
  critical_points[simplex_id] = cp;
}

template<typename T>
std::unordered_map<int, critical_point_t>
compute_critical_points(const T * U, const T * V, int r1, int r2){
  std::cout << r1 << " " << r2 << std::endl;
  // check cp for all cells
  int indices[3] = {0};
  double X1[3][2] = {
    {0, 0},
    {0, 1},
    {1, 1}
  };
  double X2[3][2] = {
    {0, 0},
    {1, 0},
    {1, 1}
  };
  double v[3][2] = {0};
  std::unordered_map<int, critical_point_t> critical_points;
  for(int i=1; i<r1-2; i++){
    if(i%100==0) std::cout << i << " / " << r1-1 << std::endl;
    for(int j=1; j<r2-2; j++){
      ptrdiff_t cell_offset = 2*(i * (r2-1) + j);
      // ptrdiff_t cell_offset = 2*(i * r2 + j);
      indices[0] = i*r2 + j;
      indices[1] = (i+1)*r2 + j;
      indices[2] = (i+1)*r2 + (j+1); 
      // cell index 0
      for(int p=0; p<3; p++){
        v[p][0] = U[indices[p]];
        v[p][1] = V[indices[p]];
      }
      check_simplex_seq(v, X1, indices, i, j, cell_offset, critical_points);
      // cell index 1
      indices[1] = i*r2 + (j+1);
      v[1][0] = U[indices[1]], v[1][1] = V[indices[1]];     
      check_simplex_seq(v, X2, indices, i, j, cell_offset + 1, critical_points);
    }
  }
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
  int DW = atoi(argv[3]);
  int DH = atoi(argv[4]);
  int sos = 1;

  std::string cp_prefix = "origin_";
  bool cp_file = false;// file_exists(cp_prefix + "_sid.dat");
  cp_file ? printf("Critical point file found!\n") : printf("Critical point Not found, recomputing\n");
  auto critical_points_0 = cp_file ? read_criticalpoints(cp_prefix) : compute_critical_points(u, v, DH, DW);

  free(u);
  free(v);

  std::cout << "read decompressed data\n";
  std::string fn_u = std::string(argv[1]) + ".out";
  std::string fn_v = std::string(argv[2]) + ".out";
  u = readfile<float>(fn_u.c_str(), num);
  v = readfile<float>(fn_v.c_str(), num);
  std::cout << "read " << num << " decompressed data done\n";

  struct timespec start, end;
  int err = 0;
  err = clock_gettime(CLOCK_REALTIME, &start);
  auto critical_points_1 = compute_critical_points(u, v, DH, DW);
  err = clock_gettime(CLOCK_REALTIME, &end);
  std::cout << "extraction time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << std::endl;
  fprintf(stderr, "#critical_points = %lu\n", critical_points_0.size());
  fprintf(stderr, "#decompressed_critical_points = %lu\n", critical_points_1.size());

  free(u);
  free(v);

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
        std::cout << "Simplex " << cp.simplex_id << " Change from " << get_critical_point_type_string(cp.type)
          << " to " << get_critical_point_type_string(cp_1.type) <<
           ", positions from " << cp.x[0] << ", " << cp.x[1] << " to " << cp_1.x[0] << ", " << cp_1.x[1] << std::endl;
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
  std::cout << "********* "<< std::endl;

  // record_criticalpoints(cp_prefix, origin, true);  

  {
    std::cout << "cp file: " << cp_file << std::endl;
    if(!cp_file) record_criticalpoints(std::string("origin_M"), origin, true);
    std::string prefix = std::string(argv[5]);
    if(m.size()) record_criticalpoints(prefix+"_M", m);    
    if(fp.size()) record_criticalpoints(prefix+"_FP", fp);    
    if(fn.size()) record_criticalpoints(prefix+"_FN", fn);    
    if(ft.size()) record_criticalpoints(prefix+"_FT", ft);    
  }
  return 0;
}