/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet
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
#ifndef TRIQS_GF_MATSUBARA_FREQ_H
#define TRIQS_GF_MATSUBARA_FREQ_H
#include "./tools.hpp"
#include "./gf.hpp"
#include "./local/tail.hpp"
#include "./domains/matsubara.hpp"
#include "./meshes/linear.hpp"
namespace triqs { namespace gfs {

 struct imfreq {};

 template<typename Opt> struct gf_mesh<imfreq,Opt> : linear_mesh<matsubara_domain<true>> {
  typedef  linear_mesh<matsubara_domain<true>> B;
  static double m1(double beta) { return std::acos(-1)/beta;}
  gf_mesh() = default;
  gf_mesh (typename B::domain_t const & d, int Nmax = 1025) : 
   B(d, d.statistic==Fermion?m1(d.beta):0, d.statistic==Fermion?(2*Nmax+1)*m1(d.beta): 2*Nmax*m1(d.beta), Nmax, without_last){}
  // use delegating...
  gf_mesh (double beta, statistic_enum S, int Nmax = 1025) :
   //gf_mesh({beta,S}, Nmax){} 
   B(typename B::domain_t(beta,S), S==Fermion?m1(beta):0, S==Fermion?(2*Nmax+1)*m1(beta): 2*Nmax*m1(beta), Nmax, without_last){}
 };

 namespace gfs_implementation { 

  //singularity
  template<typename Opt> struct singularity<imfreq,matrix_valued,Opt>  { typedef local::tail type;};
  template<typename Opt> struct singularity<imfreq,scalar_valued,Opt>  { typedef local::tail type;};

  //h5 name
  template<typename Opt> struct h5_name<imfreq,matrix_valued,Opt>      { static std::string invoke(){ return "ImFreq";}};

  /// ---------------------------  evaluator ---------------------------------
  template<typename Opt, typename Target>
   struct evaluator<imfreq,Target,Opt> {
    static constexpr int arity = 1;
     evaluator() = default;
     evaluator(int n1, int n2){} 
    typedef typename std::conditional < std::is_same<Target, matrix_valued>::value, arrays::matrix_view<std::complex<double>>, std::complex<double>>::type rtype; 
    template<typename G>
     rtype operator() (G const * g, long n)  const {return g->data()(n, arrays::ellipsis()); }
    // crucial because the mesh_point is cast in a complex, not an int !
    template<typename G>
     rtype operator() (G const * g, linear_mesh<matsubara_domain<true>>::mesh_point_t const & p)  const { return (*this)(g,p.index());}
    template<typename G>
     local::tail_view operator()(G const * g, freq_infty const &) const {return g->singularity();}
   };

  /// ---------------------------  data access  ---------------------------------
  template<typename Opt> struct data_proxy<imfreq,matrix_valued,Opt> : data_proxy_array<std::complex<double>,3> {};
  template<typename Opt> struct data_proxy<imfreq,scalar_valued,Opt> : data_proxy_array<std::complex<double>,1> {};

  // -------------------------------   Factories  --------------------------------------------------

  template<typename Opt> struct factories<imfreq,matrix_valued,Opt> : factories_one_var<imfreq,matrix_valued,Opt> {};
  template<typename Opt> struct factories<imfreq,scalar_valued,Opt> : factories_one_var<imfreq,scalar_valued,Opt> {};
 } // gfs_implementation 

}}
#endif
