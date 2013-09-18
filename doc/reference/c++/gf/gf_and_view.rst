.. highlight:: c

.. _gf_and_view:

gf and gf_view 
=================

gf is the Green function container, with its view gf_view,  defined as ::
  
  template<typename Variable, typename Target=matrix_valued, typename Opt=void> class gf;      // regular type
  template<typename Variable, typename Target=matrix_valued, typename Opt=void> class gf_view; // view type   

In this section, we describe all common properties of gf. 
In practice, one works with one of the specialization described below (:ref:`gf_special`), 
with their specific aspects (domains, interpolation methods, ....).
 

Template parameters
----------------------------

* **Variable** determines the domain of the Green function
* **Target** determines the value of the Green function
* **Opt** is currently not used [default to void]. It may be used to differentiate different Green functions
   with the same Variable, Target but different behaviour (e.g. interpolation, ...).

They can have the following values : 

+--------------------------+--------------------------------------------+
| Variable                 | Meaning                                    |
+==========================+============================================+
| gf_imfreq                | Imaginary Matsubara frequency              |
+--------------------------+--------------------------------------------+
| imtime                   | Imaginary Matsubara time                   |
+--------------------------+--------------------------------------------+
| refreq                   | Real frequency                             |
+--------------------------+--------------------------------------------+
| retime                   | Real time                                  |
+--------------------------+--------------------------------------------+
| legendre                 | Legendre polynomial representation         |
+--------------------------+--------------------------------------------+
| block_index              | Block Green function                       |
+--------------------------+--------------------------------------------+
| cartesian_product<Gs...> | Cartesian product of gf<Gs> ... functions. |
+--------------------------+--------------------------------------------+

and 

+-------------------------+-----------------------------------------------------+
| Target                  | Meaning                                             |
+=========================+=====================================================+
| matrix_valued [default] | The function is matrix valued.                      |
+-------------------------+-----------------------------------------------------+
| scalar_valued           | The function is scalar valued (double, complex...). |
+-------------------------+-----------------------------------------------------+

.. _gf_special:

Specializations
-------------------

+-----------------+------------------+--------------------+----------------------------+
| Variable/Target | matrix_valued    | scalar_valued      | G (another Green function) |
+=================+==================+====================+============================+
| imfreq          | :doc:`gf_imfreq` | :doc:`gf_imfreq_s` |                            |
+-----------------+------------------+--------------------+----------------------------+
| imtime          | :doc:`gf_imtime` | :doc:`gf_imtime_s` |                            |
+-----------------+------------------+--------------------+----------------------------+
| refreq          | :doc:`gf_refreq` | :doc:`gf_refreq_s` |                            |
+-----------------+------------------+--------------------+----------------------------+
| retime          | :doc:`gf_retime` | :doc:`gf_retime_s` |                            |
+-----------------+------------------+--------------------+----------------------------+
| block_index     |                  |                    | :doc:`block_gf<gf_block>`  |
+-----------------+------------------+--------------------+----------------------------+

.. toctree::
   :hidden:
   :maxdepth: 1

   gf_imfreq
   gf_imfreq_s
   gf_refreq
   gf_imtime
   gf_retime
   gf_block
   gf_prod


Member types 
--------------------------------------

+----------------+-------------------------------------------------------------+
| Member type    | Definitions                                                 |
+================+=============================================================+
| view_type      | The corresponding view type                                 |
+----------------+-------------------------------------------------------------+
| regular_type   | The corresponding regular type i.e. the container itself    |
+----------------+-------------------------------------------------------------+
| mesh_t         | The mesh                                                    |
+----------------+-------------------------------------------------------------+
| target_shape_t | Type storing the information to construct the target space, |
|                | Depends on the specialisation (a shape for matrix_valued    |
|                | gf, empty for a scalar valued, ... Cf Specialisations)      |
+----------------+-------------------------------------------------------------+
| data_t         | Type of the data array                                      |
+----------------+-------------------------------------------------------------+
| singularity_t  | Type of the singularity (tail, nothing...)                  |
+----------------+-------------------------------------------------------------+
| symmetry_t     | Symmetry (unused at this stage).                            |
+----------------+-------------------------------------------------------------+

Member functions
---------------------

+-------------------------------------------+------------------------------------------+
| Member function                           | Meaning                                  |
+===========================================+==========================================+
| :ref:`(constructor)<gf_constr>`           |                                          |
+-------------------------------------------+------------------------------------------+
| (destructor)                              |                                          |
+-------------------------------------------+------------------------------------------+
| :ref:`operator ()<gf_call>`               | Evaluation on a point of the domain      |
+-------------------------------------------+------------------------------------------+
| :ref:`operator []<gf_subscript>`          | Access to the value on the mesh          |
+-------------------------------------------+------------------------------------------+
| :ref:`mesh<gf_data>`                      | Access to the mesh                       |
+-------------------------------------------+------------------------------------------+
| :ref:`singularity<gf_data>`               | Access to the singularity                |
+-------------------------------------------+------------------------------------------+
| :ref:`symmetry<gf_data>`                  | Access to the symmetry                   |
+-------------------------------------------+------------------------------------------+
| :ref:`data<gf_data>`                      | Direct view of the data [Advanced]       |
+-------------------------------------------+------------------------------------------+
| :ref:`operator =<gf_reg_assign>`          | Assigns values to the container          |
+-------------------------------------------+------------------------------------------+
| :ref:`operator +=,-=,*=,/=<gf_comp_ops>`  | compound assignment operators            |
+-------------------------------------------+------------------------------------------+


.. toctree::

  :hidden:

   gf_constructors
   gf_data
   gf_assign
   gf_call
   gf_subcript
   compound_ops
   call 
   resize
   STL 

Non-member functions
------------------------


+---------------------------------+-------------------------------------------+
| Member function                 | Meaning                                   |
+=================================+===========================================+
| :ref:`swap<gf_swap>`            | Swap of 2 containers                      |
+---------------------------------+-------------------------------------------+
| :ref:`operator\<\<<gf_stream>`  | Writing to stream                         |
+---------------------------------+-------------------------------------------+


.. toctree::
   :hidden:

   stream
   swap
