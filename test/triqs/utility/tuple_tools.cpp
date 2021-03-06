/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include <triqs/utility/tuple_tools.hpp>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <functional>
#include <triqs/utility/mini_vector.hpp>

struct fun { 
 double operator()(int i, double x, double y, int k) {  return 6*k + i - 1.3*x + 2*y;}
};

struct fun2 { 
 double operator()(double x, double y) {  return x+y;}
};

struct print_t {
  template<typename T>
  void operator()(T x) { std::cerr << x << " "; }
  std::string s;
};

struct A {
template<typename T1, typename T2>  
  std::string str(T1 const & x, T2 const &y) const { std::stringstream fs; fs << "the string is "<< x<<" " << y; return fs.str();}
};

std::string my_print_str(int x, int y) { std::stringstream fs; fs << "the string is "<< x<< " " << y; return fs.str();}

namespace triqs { namespace tuple { 

}}

int main(int argc, char **argv) {

 auto t = std::make_tuple(1,2.3,4.3,8);
 auto t2 = std::make_tuple(1,2,3,4);
 auto t1 = std::make_tuple(1,2.3,4.3,8);

 {
  triqs::tuple::for_each(t, print_t());
  std::cerr << std::endl;
 }

 {
  auto res = triqs::tuple::apply(fun(),t);
  std::cerr  << " f(t) =" << res << std::endl ;
  if ( std::abs((res -  fun()(1,2.3,4.3,8))) > 1.e-13) throw std::runtime_error(" ");
 }

 {
  auto r = triqs::tuple::apply_on_zip(fun2(),t1,t2);
  std::cerr  << " [f(a,b) for (a,b) in zip(t1,t2)] =" 
   << std::get<0>(r) << " "
   << std::get<1>(r) << " "
   << std::get<2>(r) << " "
   << std::get<3>(r) << std::endl;
 }

 {
  auto res = triqs::tuple::fold([](double x,double r) { return x+r;}, t, 0);
  std::cerr  << " " << res << std::endl ;
  if ( std::abs((res -  15.6)) > 1.e-13) throw std::runtime_error(" ");
 }

 {
  auto res = triqs::tuple::fold_on_zip([](double x, double y, double r) { return x+ 2*y +r;}, t1,t2, 0);
  std::cerr  << " " << res << std::endl ;
  if ( std::abs((res -  35.6)) > 1.e-13) throw std::runtime_error(" ");
 }

 {
  auto res = triqs::tuple::apply(my_print_str, std::make_tuple(1,2));
  std::cerr  << " " << res << std::endl ;
 }

 {
  using namespace std::placeholders;
  A a;
  auto res = triqs::tuple::apply(std::bind(&A::str<int,int>,a,_1,_2), std::make_tuple(1,2));
  res = triqs::tuple::apply(std::mem_fn(&A::str<int,int>), std::make_tuple(a,1,2));
  std::cerr  << " " << res << std::endl ;
 }

 { // to mini_vector

  auto t = std::make_tuple(1,2,3.4);
  auto m = triqs::utility::tuple_to_mini_vector<double>(t);
  std::cout  << m<< std::endl ;

 }

#ifdef TRIQS_COMPILER_IS_C11_COMPLIANT

 { // filter
  std::cout  << "  ----- filter ----"<< std::endl ;
  auto t= std::make_tuple(0,1,2,3,4,"=5");
  std::cout << "filter "<< t << triqs::tuple::filter<0,2,3>(t)<< std::endl; 
  std::cout << "filter "<< t << triqs::tuple::filter<1,3,5>(t)<< std::endl; 

  std::cout << "filter out "<< t << triqs::tuple::filter_out<0,2,3>(t)<< std::endl; 
  std::cout << "filter out "<< t << triqs::tuple::filter_out<1,3,5>(t)<< std::endl; 

  auto t2= std::make_tuple(0,1);
  std::cout << "filter out "<< t2 << triqs::tuple::filter_out<0>(t2)<< std::endl; 

  typedef typename triqs::tuple::filter_t_tr< decltype(t), 0,2,3>::type TY;
  static_assert(std::is_same<TY, decltype(triqs::tuple::filter<0,2,3>(t))>::value, "EEE");
 }

 { // filter
  std::cout  << "  ----- inverse filter ----"<< std::endl ;
  auto t= std::make_tuple(1,4,5);
  auto s = std::string{"--"};
  {
   auto r = triqs::tuple::inverse_filter<6,0,2,3>(t, s);
   std::cout << "inverse filter "<< t << r<< triqs::tuple::filter<0,2,3>(r)<<std::endl; 
  }
  { auto r = triqs::tuple::inverse_filter<6,0,2,5>(t, s);
   std::cout << "inverse filter "<< t << r<< triqs::tuple::filter<0,2,5>(r)<<std::endl; 
  }
  { auto r = triqs::tuple::inverse_filter<8,0,2,5>(t, s);
   std::cout << "inverse filter "<< t << r<< triqs::tuple::filter<0,2,5>(r)<<std::endl; 
  }
  {
   auto r = triqs::tuple::inverse_filter_out<0,2,3>(t, s);
   std::cout << "inverse filter out "<< t << r<< triqs::tuple::filter_out<0,2,3>(r)<<std::endl; 
  }
  { auto r = triqs::tuple::inverse_filter_out<0,2,3,5>(t, s);
   std::cout << "inverse filter out "<< t << r<< triqs::tuple::filter_out<0,2,3,5>(r)<<std::endl; 
  }
  { auto r = triqs::tuple::inverse_filter_out<0,2,3,6>(t, s);
   std::cout << "inverse filter out "<< t << r<< triqs::tuple::filter_out<0,2,3,6>(r)<<std::endl; 
  }
 }


 { // replace
  std::cout  << "  ----- filter ----"<< std::endl ;
  auto t= std::make_tuple(0,1,2,3,4,"=5");
  auto s = std::string{"--"};
  std::cout << "replace 0,2,3"<< t << triqs::tuple::replace<0,2,3>(t,s)<< std::endl; 
  std::cout << "replace 1,3,5"<< t << triqs::tuple::replace<1,3,5>(t,s)<< std::endl; 
 }

#endif

}

