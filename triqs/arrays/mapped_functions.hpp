/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by O. Parcollet
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
#ifndef TRIQS_ARRAYS_MAPPED_FNT_H
#define TRIQS_ARRAYS_MAPPED_FNT_H
#include "./functional/map.hpp"
#include <boost/preprocessor/seq/for_each.hpp>
namespace triqs { namespace arrays {

 //C++14 will simply be ... 
 //template <typename A> decltype(auto) abs(A && a) { return map( [](auto const &x) { using std::abs; return abs(a);}, std::forward<A>(a));}

 using std::pow;
 template<typename T> TYPE_ENABLE_IF(T, std::is_integral<T>) pow(T x, int n) { return (n==0 ? 1 : pow(x,n-1)*x);}
 struct pow_wrap { 
   int n; 
   template<typename A> A operator()( A const & a) const { return pow(a,n);}
 };

 template <typename A>
  typename boost::lazy_enable_if_c<!ImmutableMatrix<A>::value,std::result_of<map_impl<pow_wrap,1>(A)>>::type
  pow(A && a, int n) { return map(pow_wrap{n})(std::forward<A>(a)); }

#define MAP_IT(FNT) \
 using std::FNT;\
 struct __triqs_##FNT##_wrap { template<typename A> auto operator()(A const & a) const DECL_AND_RETURN(FNT(a)); };\
 template <typename A> \
 typename boost::lazy_enable_if_c<ImmutableCuboidArray<A>::value,std::result_of<map_impl<__triqs_##FNT##_wrap,1>(A)>>::type\
 FNT(A && a) { return map(__triqs_##FNT##_wrap{})(std::forward<A>(a)); }
 
#define TRIQS_ARRAYS_MATH_FNT (abs)(real)(imag)(floor)

#define AUX(r, data, elem) MAP_IT(elem)
 BOOST_PP_SEQ_FOR_EACH(AUX , nil , TRIQS_ARRAYS_MATH_FNT);
#undef AUX  
#undef MAP_IT
#undef TRIQS_ARRAYS_MATH_FNT

// Functions only defined for vector, array but NOT for matrix. They act element wise.
#define MAP_IT(FNT) \
 using std::FNT;\
 struct __triqs_##FNT##_wrap { template<typename A> auto operator()(A const & a) const DECL_AND_RETURN(FNT(a)); };\
 template <typename A> \
 typename boost::lazy_enable_if_c<ImmutableArray<A>::value||ImmutableVector<A>::value,std::result_of<map_impl<__triqs_##FNT##_wrap,1>(A)>>::type\
 FNT(A && a) { return map(__triqs_##FNT##_wrap{})(std::forward<A>(a)); }
 
#define TRIQS_ARRAYS_MATH_FNT (exp)(cos)(sin)(tan)(cosh)(sinh)(tanh)(acos)(asin)(atan)(log)(sqrt)

#define AUX(r, data, elem) MAP_IT(elem)
 BOOST_PP_SEQ_FOR_EACH(AUX , nil , TRIQS_ARRAYS_MATH_FNT);
#undef AUX  
#undef MAP_IT
#undef TRIQS_ARRAYS_MATH_FNT

}}//namespace triqs::arrays 
#endif

