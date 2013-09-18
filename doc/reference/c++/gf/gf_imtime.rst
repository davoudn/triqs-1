.. highlight:: c

.. _gf_imtime: 

gf<imtime> 
===================================================

This is a specialisation of :ref:`gf_and_view` for imaginary Matsubara time.
 
Domain & mesh
----------------


Singularity
-------------

Factories
-------------

Code ::
   
   
  make_gf(MeshType && m, tqa::mini_vector<size_t,2> shape, local::tail_view const & t)
  make_gf(double beta, statistic_enum S,  tqa::mini_vector<size_t,2> shape, size_t Nmax=1025, mesh_kind mk= half_bins)
  make_gf(double beta, statistic_enum S, tqa::mini_vector<size_t,2> shape, size_t Nmax, mesh_kind mk, local::tail_view const & t)


Interpolation method
---------------------

Data storage
---------------


HDF5 storage convention
---------------------------



Examples
---------

.. compileblock::

    #include <triqs/gfs.hpp>
    using namespace triqs::gfs; using triqs::clef::placeholder;
    int main(){
     double beta=10, a = 1;
     int n_times=5;

     // First give information to build the mesh, second to build the target
     auto GF1  = gf<imtime> { {beta,Fermion,n_times}, {1,1} };
     // or a more verbose/explicit form ...
     auto GF2  = gf<imtime> { gf_mesh<imtime>{beta,Fermion,n_times}, make_shape(1,1) };
 
     // Filling the gf with something...
     placeholder<0> tau_;
     GF1(tau_) << exp ( - a * tau_) / (1 + exp(- beta * a));

     // evaluation at n=3
     std::cout << GF1(3) << " == "<<   exp ( - a * 3) / (1 + exp(- beta * a)) << std::endl;
     // the high frequency expansion was automatically computed.
     std::cout << GF1.singularity() << std::endl;
    }




