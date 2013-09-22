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
#ifndef TRIQS_GF_EVALUATOR_H
#define TRIQS_GF_EVALUATOR_H
#include "./tools.hpp"
#include "./gf.hpp"

namespace triqs { namespace gfs { 
 namespace gfs_implementation { 

 // simple evaluation : take the point on the grid...
 struct evaluator_grid_simple { 
  long n; 
  evaluator_grid_simple() = default;

  template<typename MeshType, typename PointType> 
   evaluator_grid_simple (MeshType const & m, PointType const & p) { n=p; }
  template<typename F> auto operator()(F const & f) const DECL_AND_RETURN(f (n));
 };

 // a linear interpolation 
 struct evaluator_grid_linear_interpolation { 
  double w1, w2; size_t n1, n2; 

  evaluator_grid_linear_interpolation() = default;

  template<typename MeshType, typename PointType> 
   evaluator_grid_linear_interpolation (MeshType const & m, PointType const & p, double prefactor=1) { // delegate !
    bool in; double w;
    std::tie(in, n1, w) = windowing(m,p);
    //std::cout  << in << " "<< n1 << " "<< w << " " << p << std::endl;
    if (!in) TRIQS_RUNTIME_ERROR <<" Evaluation out of bounds";
    w1 = prefactor * (1-w); w2 = prefactor * w; n2 = n1 +1;
   }

  template<typename F> auto operator()(F const & f) const DECL_AND_RETURN(w1 * f(n1) + w2 * f (n2));
 };

 // the evaluator for various types.
 template<typename MeshType> struct evaluator_fnt_on_mesh;

 // can not use inherited constructors, too recent...
#define TRIQS_INHERIT_AND_FORWARD_CONSTRUCTOR(NEWCLASS,CLASS) : CLASS { template<typename ...T> NEWCLASS(T &&... t) : CLASS(std::forward<T>(t)...){};}; 

 //template<> struct evaluator_fnt_on_mesh<imfreq> TRIQS_INHERIT_AND_FORWARD_CONSTRUCTOR(evaluator_fnt_on_mesh, evaluator_grid_simple);
 //template<> struct evaluator_fnt_on_mesh<imtime> TRIQS_INHERIT_AND_FORWARD_CONSTRUCTOR(evaluator_fnt_on_mesh, evaluator_grid_linear_interpolation);
 //template<> struct evaluator_fnt_on_mesh<retime> TRIQS_INHERIT_AND_FORWARD_CONSTRUCTOR(evaluator_fnt_on_mesh, evaluator_grid_linear_interpolation);
 //template<> struct evaluator_fnt_on_mesh<refreq> TRIQS_INHERIT_AND_FORWARD_CONSTRUCTOR(evaluator_fnt_on_mesh, evaluator_grid_linear_interpolation);

 // 
 template<typename Variable>
  struct evaluator_one_var { 
   public :
    static constexpr int arity = 1;
    evaluator_one_var() = default;
    evaluator_one_var(int n1, int n2) {}; // SUPRESS THIS : _tmp(n1,n2) {} // WHAT happen in resize ??
#ifndef TRIQS_COMPILER_OBSOLETE_GCC    
 template<typename G> auto operator()(G const * g, double x) const DECL_AND_RETURN(evaluator_fnt_on_mesh<Variable> (g->mesh(),x)(on_mesh(*g)));
#else
 template<typename G> auto operator()(G const * g, double x) const DECL_AND_RETURN(evaluator_fnt_on_mesh<Variable> (g->mesh(),x)(typename G::_on_mesh_wrapper_const(*g)));
#endif

    template<typename G> typename G::singularity_t const & operator()(G const * g,freq_infty const &) const {return g->singularity();}
  };

 }
}}
#endif
