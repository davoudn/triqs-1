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
#ifndef TRIQS_GF_BLOCK_H
#define TRIQS_GF_BLOCK_H
#include "./tools.hpp"
#include "./gf.hpp"
#include "./local/tail.hpp"
#include "./meshes/discrete.hpp"

namespace triqs { namespace gfs {

 struct block_index {};

 template<typename Opt> struct gf_mesh<block_index,Opt> : discrete_mesh<discrete_domain> {
  typedef discrete_mesh<discrete_domain> B;
  gf_mesh() = default;
  gf_mesh(size_t s) : B(s) {}
  gf_mesh(discrete_domain const & d) : B(d) {}
  gf_mesh(std::initializer_list<std::string> const & s) : B(s){}
 };

 namespace gfs_implementation { 

  template<typename Target, typename Opt> struct h5_name<block_index,Target,Opt>      { static std::string invoke(){ return  "BlockGf";}};

  /// ---------------------------  h5_rw ---------------------------------

  template <typename Target,typename Opt> struct h5_ops<block_index,Target,Opt> { 

   template<typename DataType, typename GF> 
    static void write(h5::group g, std::string const & s, DataType const & data, GF const & gf) {
     auto gr =  g.create_group(s);
     for (size_t i =0; i<gf.mesh().size(); ++i) h5_write(gr,gf.mesh().domain().names()[i],data[i]);
    }

   template<typename DataType,typename GF> 
    static void read(h5::group g, std::string const & s, DataType & data, GF const & gf) {
     auto gr =  g.create_group(s);
     for (size_t i =0; i<gf.mesh().size(); ++i) h5_read(gr,gf.mesh().domain().names()[i],data[i]);
    }
  };

  /// ---------------------------  data access  ---------------------------------

  template<typename Target, typename Opt> struct data_proxy<block_index,Target,Opt> : data_proxy_vector <typename regular_type_if_exists_else_type<Target>::type>{};

  // -------------------------------   Factories  --------------------------------------------------

  template<typename Target, typename Opt>
   struct factories<block_index,Target,Opt>  { 
    typedef gf_mesh<block_index, Opt> mesh_t;
    typedef gf<block_index,Target> gf_t;
    typedef gf_view<block_index,Target> gf_view_t;
    typedef std::initializer_list<Target> target_shape_t;

    // target_shape here is a init list of gf : if {G} -> repeat the same GF
    static typename gf_t::data_t make_data(mesh_t const & m, target_shape_t const & t_shape) { 
     if ((t_shape.size() !=1) && (t_shape.size()!=m.size())) TRIQS_RUNTIME_ERROR<<" gf block construct size mismatch";
     std::vector<Target> r; 
     for (auto const & g : t_shape) r.push_back(g);
     for (int i= t_shape.size(); i<m.size(); ++i) r.push_back(r[0]);
     return r;
    }

    static typename gf_t::singularity_t       make_singularity (mesh_t const & m, target_shape_t) { return {};}
    static evaluator<block_index, Target,Opt> make_evaluator   (mesh_t const & m, target_shape_t) { return {};}
   };

 } // gfs_implementation

 // -------------------------------   Free functions   --------------------------------------------------

 template<typename G0, typename ... G> 
   gf_view<block_index, typename std::remove_reference<G0>::type::view_type> make_block_gf_view(G0 && g0, G && ... g) { 
    auto V = std::vector<typename std::remove_reference<G0>::type::view_type>{std::forward<G0>(g0), std::forward<G>(g)...};
    return { {V.size()}, std::move(V), nothing(), nothing() } ;
    //return { gf_mesh<block_index, Opt> {V.size()}, std::move(V), nothing(), nothing() } ;
   }

 template<typename GF>
  static gf<block_index,GF> make_block_gf_from_vector (std::vector<GF> V)      { return { {V.size()}, std::move(V), nothing(), nothing()} ; }

 template<typename GF>
  static gf<block_index,GF> make_block_gf_from_name_and_vector (std::vector<std::string> const & block_names, std::vector<GF> V)  { 
   return { {block_names}, std::move(V), nothing(), nothing()}; 
  }

 template<typename GF, typename GF2>
  static gf_view<block_index,GF> make_block_gf_view_from_vector (std::vector<GF2> V)      { return { {V.size()}, std::move(V), nothing(), nothing()} ; }

 // a simple function to get the number of blocks
 template<typename T> size_t n_blocks (gf<block_index,T> const & g)      { return g.mesh().size();}
 template<typename T> size_t n_blocks (gf_view<block_index,T> const & g) { return g.mesh().size();}

#ifndef TRIQS_COMPILER_IS_OBSOLETE
 template<typename T> using block_gf = gf<block_index, gf<T>>;
#endif

 // also experimental
 // an iterator over the block
 template<typename Target, typename Opt>
  class block_gf_iterator : 
   public boost::iterator_facade< block_gf_iterator<Target,Opt>, typename Target::view_type , boost::forward_traversal_tag, typename Target::view_type  > {
    friend class boost::iterator_core_access;
    typedef gf_view<block_index,Target,Opt> big_gf_t;
    typedef typename big_gf_t::mesh_t::const_iterator mesh_iterator_t;
    big_gf_t big_gf;
    mesh_iterator_t mesh_it;

    typename Target::view_type const & dereference() const { return big_gf[*mesh_it];}
    bool equal(block_gf_iterator const & other) const { return ((mesh_it == other.mesh_it));}
    public:
    block_gf_iterator(gf_view<block_index,Target,Opt> bgf, bool at_end = false): big_gf(std::move(bgf)), mesh_it(&big_gf.mesh(),at_end) {}
    void increment() { ++mesh_it; }
    bool at_end() const { return mesh_it.at_end();}
   };

 template<typename Target, typename Opt, bool B>
  block_gf_iterator<Target,Opt> begin(gf_impl<block_index,Target,Opt,B> const & bgf) { return {bgf,false};}

 template<typename Target, typename Opt, bool B>
  block_gf_iterator<Target,Opt> end(gf_impl<block_index,Target,Opt,B> const & bgf) { return {bgf,true};}

}}
#endif


