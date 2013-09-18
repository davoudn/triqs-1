.. highlight:: c

.. _gf_constr:

Constructors
====================

Constructors of gf
---------------------------


+--------------------------------+-------------------------------------------------------------------------------------------------------------------+
| Constructors                   | Comments                                                                                                          |
+================================+===================================================================================================================+
| gf()                           | Empty gf                                                                                                          |
+--------------------------------+-------------------------------------------------------------------------------------------------------------------+
| gf(const gf &)                 | Copy construction                                                                                                 |
+--------------------------------+-------------------------------------------------------------------------------------------------------------------+
| gf(gf &&)                      | Move construction                                                                                                 |
+--------------------------------+-------------------------------------------------------------------------------------------------------------------+
| gf(gf_view & g)                | Make a clone of the view.                                                                                         |
+--------------------------------+-------------------------------------------------------------------------------------------------------------------+
| gf(mesh_t m, target_shape_t s) | Constructs a gf with mesh m and a target shape s. s depends on the specialization: it is                          |
|                                | the information needed to construct the gf, aside from the mesh itself : typically the                            |
|                                | size of the matrix for a matrix_valued gf. target_shape_t is described in the                                     |
|                                | specialization section. section XXX). The function is initialized to 0                                            |
+--------------------------------+-------------------------------------------------------------------------------------------------------------------+
| gf(mesh_t m,                   | *[Advanced]* Construct a gf from its elements (a mesh, a data array, a singularity, a symmetry, and an evaluator. |
| data_t dat,                    | Normally should be reserved for lib purposes...(remove from the doc ?).                                           |
| singularity_view_t const & si, |                                                                                                                   |
| symmetry_t const & s,          |                                                                                                                   |
| evaluator_t const & eval       |                                                                                                                   |
+--------------------------------+-------------------------------------------------------------------------------------------------------------------+



Examples
------------

.. compileblock::

    #include <triqs/gfs.hpp>
    using triqs::gfs::gf; using triqs::gfs::matrix; 
    using triqs::gfs::vector; using triqs::gfs::permutation; 
    int main(){
      
      // A 3d gf of long, C ordering, no option
      gf<long, 3> A3(1,2,3);
      
      // A 2d gf of double, C ordering, with explicit Bound Checking
      gf<double, 2> B(1,2);

      // a matrix of long
      matrix<long> M(2,2);
      
      // a vector of double
      vector<double> V(10);

      // gfs with custom TraversalOrder  

      // C-style
      gf<long, 3, 0, permutation(2,1,0)> A0(2,3,4);       
      gf<long, 3, 0> A0b; // same type but empty      
     
      // Fortran-style
      gf<long, 3, TRAVERSAL_ORDER_FORTRAN> A4 (2,3,4);
      gf<long, 3, 0, permutation(0,1,2)> A1b; //same type but empty      

      // custom :  (i,j,k)  : index j is fastest, then k, then i
      gf<long, 3, 0, permutation(1,0,2)> A2(2,3,4); 
    }
   

-------------------------------

gf/gf_view have very basic constructors : 
default, copy, move, and one constructor from the data used by the functions (reserved for advanced users).

Various specializations however provides several factories, adapted to each case, of the form ::

  auto g= make_gf<Variable, Target, Opt> ( ....) ;

This is the recommended way to construct `gf` objects.
Cf examples in various specializations.


