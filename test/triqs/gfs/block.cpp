#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#include <triqs/gfs.hpp> 
using namespace triqs::gfs;
#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl<<std::endl;

int main() {
 try {
  double beta =1;
  auto G1 = gf<imfreq> ({beta, Fermion}, {2,2});
  auto G2 = G1;
  auto G3 = G2;

#ifndef TRIQS_COMPILER_IS_OBSOLETE
  // construct some block functions
  auto B1 = block_gf<imfreq> (3, {G1}); 
  auto B2 = block_gf<imfreq> (3, {G1,G1,G1});
  // constructor from vector ... 
  auto B3 = block_gf<imfreq> ({"a","b","c"}, {G1,G1,G1}); 

  B1[0][0] = 98;
  //not implemented yet
  //B3["a"][0] = 98;
#endif

  auto View =  make_block_gf_view(G1,G2,G3);

  std::cout  << "Number of blocks " << View.mesh().size()<<std::endl ;
  auto g0 = View[0];
  auto g0v = View[0]();

  auto Gv = g0();

  Gv[0] = 20;
  TEST( G1( 0) ) ;
  Gv[0] = 0;

  g0v[0] = 3.2;
  TEST( G1( 0) ) ;

  // Operation
  g0[0] = 3.2;
  TEST( View[0](0) ) ;
  View = View/2; 
  TEST( View[0](0) ) ;

  // try the loop over the block.
  for (auto g : View) { g[0] = 20;}
 }
 catch(std::exception const & e) { std::cout << e.what() << std::endl;}
}
