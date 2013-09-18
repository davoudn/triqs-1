#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK

#include <triqs/gfs.hpp> 
#include <triqs/gfs/local/fourier_real.hpp> 

using namespace triqs::gfs;
using  triqs::clef::placeholder;

int main() {

 // scalar valued gf_vertex
 using gf_vertex_t        = gf<cartesian_product<imfreq,imfreq,imfreq>, scalar_valued>;
 using gf_vertex_tensor_t = gf<cartesian_product<imfreq,imfreq,imfreq>, tensor_valued<3>>;

 try { 

  double beta =10.0;
  int n_im_freq=10;

  auto m = gf_mesh<imfreq> {beta, Fermion, n_im_freq};

  auto vertex = gf_vertex_t { {m,m,m} };
 
  placeholder<0> w0_;
  placeholder<1> w1_;
  placeholder<2> w2_;
  
  vertex (w0_, w1_, w2_) << w0_ + w1_ + w2_;
  
  vertex [{0,0,0}]  = 10;
 
  auto v = on_mesh(vertex);
  v(0,0,0) *=2;

  std::cout << vertex(0,0,0)<< std::endl;
  
  //saving 
  H5::H5File file("vertex.h5", H5F_ACC_TRUNC );
  h5_write(file, "v", vertex);
  
  // now with indices
  auto vertex3 = gf_vertex_tensor_t { {m,m,m} , {2,2,2} };
  vertex3 [{0,0,0}](0,0,0)  = 10;
  std::cout << vertex3(0,0,0)<< std::endl;
  
 }
 catch(std::exception const & e ) { std::cout  << "error "<< e.what()<< std::endl;}
}
