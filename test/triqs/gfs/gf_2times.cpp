#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#include <triqs/gfs.hpp> 
#include <triqs/gfs/two_real_times.hpp> 
using namespace triqs::gfs;
using namespace triqs::arrays;
#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl<<std::endl;

int main() {

 try {

 auto G  = make_gf<two_real_times>( 10,100,make_shape(2,2));
 auto G2 = make_gf<two_real_times>( 10,100,make_shape(2,2));

 triqs::clef::placeholder<0> t_;
 triqs::clef::placeholder<1> tp_;

 array<double,2> A(2,2);
 A(t_,tp_)  << t_ - 3*tp_;
 std::cout  <<A << std::endl ;

 G(t_,tp_)  << t_ - 3*tp_;
 G2(t_,tp_) << t_ + 3*tp_;
  
 G2(t_,tp_) << 2* G(tp_,t_);

 TEST( G(1,1) );
 TEST( G[G.mesh()(1,1) ]);
 TEST( G.on_mesh(1,1));

 //G2(t_,tp_) << G(tp_,t_);
 TEST( G(2,1) );
 TEST( G2(1,2) );
 TEST( G(1,2) );
 TEST( G2(2,1) );

 //TEST( G2(2,1,3) ); // should not compile
 
 }
 catch( std::exception const &e) { std::cout  << e.what()<< std::endl;}
}
