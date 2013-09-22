//#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#include <triqs/gfs.hpp> 
using namespace triqs::gfs;
#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl<<std::endl;

int main() {
 try { 
  double beta =1;
  auto G =  gf<imfreq, scalar_valued> {{beta, Fermion}};

  std::cout << G.singularity() << std::endl ;
  triqs::clef::placeholder<0> om_;
  G(om_) << 1/(om_ + 2.3);

  // test hdf5 
  H5::H5File file("gf_scalar.h5", H5F_ACC_TRUNC);
  h5_write(file, "g", G);
  h5_write(file, "gm", reinterpret_scalar_valued_gf_as_matrix_valued(G));

 }
 catch(std::exception const & e) { std::cout << e.what() << std::endl; return 1;}


}
