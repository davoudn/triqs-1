/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012-2013 by O. Parcollet
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
#ifndef TRIQS_CLEF_CORE_H
#define TRIQS_CLEF_CORE_H
#include <triqs/utility/first_include.hpp>
#include <triqs/utility/macros.hpp>
#include <triqs/utility/compiler_details.hpp>
#include <tuple>
#include <type_traits>
#include <functional>
#include <memory>
#include <complex>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/arithmetic/inc.hpp>

#define TRIQS_CLEF_MAXNARGS 8

namespace triqs { namespace clef { 
 typedef unsigned long long ull_t;
 namespace tags { struct function_class{}; struct function{}; struct subscript{}; struct terminal{}; struct if_else{}; struct unary_op{}; struct binary_op{}; }

 // Compute the type to put in the expression tree.
 // If T is a lvalue reference, pack it into a reference_wrapper, unless force_copy_in_expr<T>::value == true
 // If T is an rvalue reference, we store it as the type (using move semantics).
 template<typename T> struct force_copy_in_expr : std::false_type{};
 template<typename T> struct force_copy_in_expr<T const> : force_copy_in_expr<T>{};

 template< class T > struct expr_storage_t      {typedef T type;};
 template< class T > struct expr_storage_t<T&>  : std::conditional<force_copy_in_expr<T>::value, T ,std::reference_wrapper<T>>{};
 template< class T > struct expr_storage_t<T&&> {typedef T type;};
 
 /* ---------------------------------------------------------------------------------------------------
  *  Placeholder and corresponding traits
  *  --------------------------------------------------------------------------------------------------- */
 template<int i, typename T> class pair; // forward

  // a placeholder is an empty struct, labelled by an int.
 template<int N> struct placeholder {
  static_assert( (N>=0) && (N<64) , "Placeholder number limited to [0,63]");
  static constexpr int index = N;
  template <typename RHS> pair<N,RHS> operator = (RHS && rhs) { return {std::forward<RHS>(rhs)};} 
 };
 
 // placeholder will always be copied (they are empty anyway).
 template< int N > struct force_copy_in_expr<placeholder<N>> : std::true_type{};

 // represent a couple (placeholder, value). 
 template<int N, typename U> struct pair {
  U rhs;
  static constexpr int p = N;
  typedef typename remove_cv_ref <U>::type value_type; 
 };

 // ph_set is a trait that given a pack of type, returns the set of placeholders they contain
 // it returns a int in binary coding : bit N in the int is 1 iif at least one T is lazy and contains placeholder<N>
 template<typename... T> struct ph_set;
 template<typename T0, typename... T> struct ph_set<T0,T...>{static constexpr ull_t value= ph_set<T0>::value| ph_set<T...>::value;};
 template<typename T> struct ph_set<T>                      {static constexpr ull_t value= 0;};
 template<int N>      struct ph_set<placeholder<N>>         {static constexpr ull_t value= 1ull<<N;};
 template<int i, typename T> struct ph_set<pair<i,T> > : ph_set<placeholder<i>>{};

 // in_any_lazy : trait to detect if any of Args is a lazy expression
 template<typename... Args>             struct is_any_lazy                     :  std::false_type {};
 template<int N>                        struct is_any_lazy <placeholder <N> >  :  std::true_type  {};
 template <typename T>                  struct is_any_lazy< T >                :  std::false_type {};
 template <typename T>                  struct is_any_lazy< T&& >              :  is_any_lazy<T> {};
 template <typename T>                  struct is_any_lazy< T& >               :  is_any_lazy<T> {};
 template <typename T>                  struct is_any_lazy< T const >          :  is_any_lazy<T> {};
 template<typename T, typename... Args> struct is_any_lazy<T, Args...>         :
  std::integral_constant<bool, is_any_lazy<T>::value || is_any_lazy<Args...>::value> {}; 

 template<typename T>
 constexpr bool ClefExpression() { return is_any_lazy<T>::value;}

 template<typename T> struct is_clef_expression : is_any_lazy<T>{};

/* ---------------------------------------------------------------------------------------------------
  * Node of the expression tree
  *  --------------------------------------------------------------------------------------------------- */
 template<typename Tag, typename... T> struct expr {
  // T can be U, U & (a reference or a value).
  typedef std::tuple<T...> childs_t; 
  childs_t childs; 
  expr(expr const & x) = default; 
  expr(expr && x) noexcept : childs(std::move(x.childs)) {}
  // a constructor with the Tag make it unambiguous with other constructors...
  template<typename... Args> expr(Tag, Args&&...args) : childs(std::forward<Args>(args)...) {} 
  // [] returns a new lazy expression, with one more layer
  template<typename Args> 
   expr<tags::subscript, expr, typename expr_storage_t<Args>::type > operator[](Args && args) const 
   { return {tags::subscript(), *this,std::forward<Args>(args)};}
  // () also ...
  template< typename... Args > 
   expr<tags::function, expr, typename expr_storage_t<Args>::type...> operator()(Args && ... args) const 
   { return {tags::function(), *this,std::forward<Args>(args)...};}
  // assignement is in general deleted
  expr & operator= (expr const &)  = delete; // no ordinary assignment
  expr & operator= (expr &&)  = default; // move assign ok 
  // however, this is ok in the case f(i,j) = expr, where f is a clef::function
  template<typename RHS, typename CH = childs_t> 
   ENABLE_IF(std::is_base_of<tags::function_class, typename std::tuple_element<0,CH>::type>)
   operator= (RHS const & rhs) { *this << rhs;}
  template<typename RHS, typename CH = childs_t> 
   DISABLE_IF(std::is_base_of<tags::function_class, typename std::tuple_element<0,CH>::type>)
   operator= (RHS const & rhs) = delete; 
 };
 // set some traits
 template<typename Tag, typename... T> struct ph_set< expr<Tag,T... > > : ph_set<T...> {};
 template<typename Tag, typename... T> struct is_any_lazy< expr<Tag,T... > >: std::true_type {};
 // if we want that subexpression are copied ?
 template<typename Tag, typename... T> struct force_copy_in_expr< expr<Tag,T... > > : std::true_type{};

 /* ---------------------------------------------------------------------------------------------------
  * The basic operations put in a template.... 
  *  --------------------------------------------------------------------------------------------------- */
 template<typename Tag> struct operation;

 // a little function to clean the reference_wrapper
 template<typename U> U & _cl(U & x) { return x;}
 template<typename U> U & _cl(std::reference_wrapper<U> x) { return x.get();}

 /// Terminal 
 template<> struct operation<tags::terminal> { template<typename L> L operator()(L&& l) const { return std::forward<L>(l);} };

 /// Function call 
 template<> struct operation<tags::function> { 
  template<typename F, typename... Args> auto operator()(F const & f, Args const & ... args) const DECL_AND_RETURN(_cl(f)(_cl(args)...));
 };

 /// [ ] Call
 template<> struct operation<tags::subscript> { 
  template<typename F, typename Args> auto operator()(F const & f, Args const & args) const DECL_AND_RETURN(_cl(f)[_cl(args)]);
 };

 // all binary operators....
#define TRIQS_CLEF_OPERATION(TAG,OP)\
 namespace tags { struct TAG : binary_op { static const char * name() { return BOOST_PP_STRINGIZE(OP);} };}\
 template<> struct operation<tags::TAG> {\
  template<typename L, typename R> auto operator()(L  const & l, R  const & r) const DECL_AND_RETURN ( _cl(l) OP _cl(r));\
 };\
 template<typename L, typename R>\
 typename std::enable_if<is_any_lazy<L,R>::value, expr<tags::TAG,typename expr_storage_t<L>::type,typename expr_storage_t<R>::type> >::type \
 operator OP (L && l, R && r) { return {tags::TAG(),std::forward<L>(l),std::forward<R>(r)};}\

 TRIQS_CLEF_OPERATION(plus,       +);
 TRIQS_CLEF_OPERATION(minus,      -);
 TRIQS_CLEF_OPERATION(multiplies, *);
 TRIQS_CLEF_OPERATION(divides,    /);
 TRIQS_CLEF_OPERATION(greater,    >);
 TRIQS_CLEF_OPERATION(less,       <);
 TRIQS_CLEF_OPERATION(leq,        <=);
 TRIQS_CLEF_OPERATION(geq,        >=);
 TRIQS_CLEF_OPERATION(eq,        ==);
#undef TRIQS_CLEF_OPERATION

 // all unary operators....
#define TRIQS_CLEF_OPERATION(TAG,OP)\
 namespace tags { struct TAG : unary_op { static const char * name() { return BOOST_PP_STRINGIZE(OP);} };}\
 template<> struct operation<tags::TAG> {\
  template<typename L> auto operator()(L const & l) const DECL_AND_RETURN (OP _cl(l));\
 };\
 template<typename L>\
 typename std::enable_if<is_any_lazy<L>::value, expr<tags::TAG,typename expr_storage_t<L>::type> >::type \
 operator OP (L && l) { return {tags::TAG(),std::forward<L>(l)};}\

 TRIQS_CLEF_OPERATION(negate,      -);
 TRIQS_CLEF_OPERATION(loginot,     !);
#undef TRIQS_CLEF_OPERATION

 /// the only ternary node :  expression if
 template<> struct operation<tags::if_else> { 
  template<typename C, typename A, typename B> auto operator()(C const & c, A const & a, B const & b) const DECL_AND_RETURN (_cl(c) ? _cl(a): _cl(b));
 };
 // operator is : if_else( Condition, A, B)
 template<typename C, typename A, typename B> 
  expr<tags::if_else,typename expr_storage_t<C>::type,typename expr_storage_t<A>::type,typename expr_storage_t<B>::type>
  if_else(C && c, A && a, B && b) { return {tags::if_else(),std::forward<C>(c),std::forward<A>(a),std::forward<B>(b)};}

 /* ---------------------------------------------------------------------------------------------------
  * Evaluation of the expression tree.
  *  --------------------------------------------------------------------------------------------------- */

 template<typename Tag, bool IsLazy,  typename... Args> struct operation2;
 template<typename Tag, typename... Args> struct operation2<Tag, true, Args...> { 
  typedef expr<Tag,typename remove_cv_ref<Args>::type ...> rtype;
  rtype operator()(Args const & ...  args) const  {return rtype {Tag(), args...};} 
 };
 template<typename Tag, typename... Args> struct operation2<Tag, false, Args...> { 
  typedef typename std::remove_reference<typename std::result_of<operation<Tag>(Args...)>::type>::type rtype;
  // remove the reference because of ternary if_else in which decltype returns a ref...
  rtype operator()(Args const & ... args) const  {return operation<Tag>()(args...); } 
 };

 // Generic case : do nothing (for the leaf of the tree except placeholder)
 template<typename T, typename ... Pairs> struct evaluator{
  typedef T rtype;
  rtype operator()(T const & k, Pairs const &... pairs) {return k;}
 };

 // placeholder
 template<int N, int i, typename T, typename... Pairs> struct evaluator< placeholder<N>, pair<i,T>, Pairs... > {
  typedef evaluator< placeholder<N>, Pairs...> eval_t;
  typedef typename eval_t::rtype rtype;
  rtype operator()(placeholder<N>, pair<i,T> const &, Pairs const& ... pairs) { return eval_t()(placeholder<N>(), pairs...);} 
 };
 template<int N, typename T, typename... Pairs> struct evaluator< placeholder<N>, pair<N,T>, Pairs... > {
  typedef T const & rtype;
  //typedef typename pair<N,T>::value_type const & rtype;
  rtype operator()(placeholder<N>, pair<N,T> const & p, Pairs const& ...) { return p.rhs;}
 };

 // general expr node
 template<typename Tag, typename... Childs, typename... Pairs> struct evaluator<expr<Tag, Childs...>, Pairs...> {
  typedef operation2<Tag, is_any_lazy<typename evaluator<Childs, Pairs...>::rtype... >::value, typename evaluator<Childs, Pairs...>::rtype... > OPTYPE;
  typedef typename OPTYPE::rtype rtype;

  // first done manually for clearer error messages ... 
  template< int arity = sizeof...(Childs)>
  typename std::enable_if< arity==1, rtype>::type
  operator()(expr<Tag, Childs...> const & ex, Pairs const & ... pairs) const
  { return OPTYPE()(eval(std::get<0>(ex.childs),pairs...) );}
 
  template< int arity = sizeof...(Childs)>
  typename std::enable_if< arity==2, rtype>::type
  operator()(expr<Tag, Childs...> const & ex, Pairs const & ... pairs) const
  { return OPTYPE()(eval(std::get<0>(ex.childs),pairs...),eval(std::get<1>(ex.childs),pairs...) );}
 
#define AUX(z,p,unused)  eval(std::get<p>(ex.childs),pairs...) 
#define IMPL(z, NN, unused)  \
  template< int arity = sizeof...(Childs)>\
  typename std::enable_if< arity==NN, rtype>::type\
  operator()(expr<Tag, Childs...> const & ex, Pairs const & ... pairs) const\
  { return OPTYPE()(BOOST_PP_ENUM(NN,AUX,nil));}
  BOOST_PP_REPEAT_FROM_TO(3,BOOST_PP_INC(TRIQS_CLEF_MAXNARGS), IMPL, nil);
#undef AUX
#undef IMPL 
 }; 

#ifdef TRIQS_CLEF_EVAL_SHORT_CIRCUIT
 // A short circuit if intersection of ph and is 0, no need to evaluate the whole tree......
 // Seems useless, && the second eval is not correct if hte expression is a terminal. 
 template<typename T, typename... Pairs> 
  typename std::enable_if< (ph_set<T>::value & ph_set<Pairs...>::value) !=0, typename evaluator<T,Pairs...>::rtype > ::type
  eval (T const & ex, Pairs const &... pairs) { return evaluator<T, Pairs...>()(ex, pairs...); }

 template<typename T, typename... Pairs> 
  typename std::enable_if< (ph_set<T>::value & ph_set<Pairs...>::value) ==0, T const &> ::type
  eval (T const & ex, Pairs const &... pairs) { return ex;}
#else
 // The general eval function for expressions
 template<typename T, typename... Pairs> 
  typename evaluator<T,Pairs...>::rtype 
  eval (T const & ex, Pairs const &... pairs) { return evaluator<T, Pairs...>()(ex, pairs...); }
#endif

 /* ---------------------------------------------------------------------------------------------------
  * make_function : transform an expression to a function  
  *  --------------------------------------------------------------------------------------------------- */

 template< typename Expr, int... Is> struct make_fun_impl { 
  Expr ex; // keep a copy of the expression
  make_fun_impl(Expr const & ex_) : ex(ex_) {} 

  // gcc 4.6 crashes (!!!) on the first variant
#ifndef TRIQS_COMPILER_OBSOLETE_GCC
  template<typename... Args> 
   typename evaluator<Expr,pair<Is,Args>...>::rtype 
   operator()(Args &&... args) const 
   { return evaluator<Expr,pair<Is,Args>...>() ( ex, pair<Is,Args>{std::forward<Args>(args)}...); }
#else
  template<typename... Args> struct __eval { 
   typedef evaluator<Expr,pair<Is,Args>...> eval_t;
   typedef typename eval_t::rtype rtype;
   rtype operator()(Expr const &ex , Args &&... args) const { return eval_t() ( ex, pair<Is,Args>{std::forward<Args>(args)}...); }
  };
  template<typename... Args> 
   typename __eval<Args...>::rtype operator()(Args &&... args) const { return __eval<Args...>() ( ex, std::forward<Args>(args)...); } 
#endif
 };

 // values of the ph, excluding the Is ...
 template<ull_t x, int... Is> struct ph_filter;
 template<ull_t x, int I0, int... Is> struct ph_filter<x,I0,Is...> { static constexpr ull_t value =  ph_filter<x,Is...>::value & (~ (1ull<<I0));};
 template<ull_t x> struct ph_filter<x> { static constexpr ull_t value = x; }; 

 template< typename Expr, int... Is> struct ph_set<make_fun_impl<Expr,Is...> >  { static constexpr ull_t value = ph_filter <ph_set<Expr>::value, Is...>::value;};
 template< typename Expr, int... Is> struct is_any_lazy<make_fun_impl<Expr,Is...> > : std::integral_constant<bool,ph_set<make_fun_impl<Expr,Is...>>::value !=0>{};
 template< typename Expr, int... Is> struct force_copy_in_expr<make_fun_impl<Expr,Is...> > : std::true_type{};

 template< typename Expr, int... Is,typename... Pairs> struct evaluator<make_fun_impl<Expr,Is...>, Pairs...> {
  typedef evaluator<Expr,Pairs...> e_t;
  typedef make_fun_impl<typename e_t::rtype, Is...> rtype;
  rtype operator()(make_fun_impl<Expr,Is...> const & f, Pairs const & ... pairs) const { return rtype( e_t()(f.ex, pairs...));}
 };

 template< typename Expr, typename ... Phs> 
  make_fun_impl<typename remove_cv_ref <Expr>::type,Phs::index...> 
  make_function(Expr && ex, Phs...) { return {ex}; }

 namespace result_of {
  template< typename Expr, typename ... Phs> struct make_function { 
   typedef make_fun_impl<typename remove_cv_ref<Expr>::type,Phs::index...> type;
  };
 }

 template<int ... N> struct ph_list {};
 template<int ... N> ph_list<N...> var( placeholder<N> ...) { return {};}

 template<typename Expr, int ... N> 
 auto operator >> (ph_list<N...>, Expr const & ex) DECL_AND_RETURN( make_function(ex, placeholder<N>()...));

 /* --------------------------------------------------------------------------------------------------
  *  make_function
  *  x_ >> expression  is the same as make_function(expression,x)
  * --------------------------------------------------------------------------------------------------- */

 template <int N, typename Expr>
  make_fun_impl<Expr,N > operator >> (placeholder<N> p, Expr&& ex) { return {ex}; } 

 /* ---------------------------------------------------------------------------------------------------
  * Auto assign for ()
  *  --------------------------------------------------------------------------------------------------- */

 // by default it is deleted = not implemented : every class has to define it...
 template<typename T, typename F> void triqs_clef_auto_assign (T,F) = delete;

 // remove the ref_wrapper, terminal ...
 template<typename T, typename F> void triqs_clef_auto_assign (std::reference_wrapper<T> R     ,F && f) { triqs_clef_auto_assign(R.get(),std::forward<F>(f));}
 template<typename T, typename F> void triqs_clef_auto_assign (expr<tags::terminal,T> const & t,F && f) { triqs_clef_auto_assign(std::get<0>(t.childs),std::forward<F>(f));}

 // auto assign of an expr ? (for chain calls) : just reuse the same operator
 template<typename Tag, typename... Childs, typename RHS> 
  void triqs_clef_auto_assign (expr<Tag,Childs...> && ex, RHS const & rhs) { ex << rhs;}

 template<typename Tag, typename... Childs, typename RHS> 
  void triqs_clef_auto_assign (expr<Tag,Childs...> const & ex, RHS const & rhs) { ex << rhs;}

 // The case A(x_,y_) = RHS : we form the function (make_function) and call auto_assign (by ADL)
 template<typename F, typename RHS, int... Is> 
  void operator<<(expr<tags::function, F, placeholder<Is>...> && ex, RHS && rhs) { 
   triqs_clef_auto_assign(std::get<0>(ex.childs), make_function(std::forward<RHS>(rhs), placeholder<Is>()...));
  }
 template<typename F, typename RHS, int... Is> 
  void operator<<(expr<tags::function, F, placeholder<Is>...> const & ex, RHS && rhs) { 
   triqs_clef_auto_assign(std::get<0>(ex.childs), make_function(std::forward<RHS>(rhs), placeholder<Is>()...));
  }
 template<typename F, typename RHS, int... Is> 
  void operator<<(expr<tags::function, F, placeholder<Is>...> & ex, RHS && rhs) { 
   triqs_clef_auto_assign(std::get<0>(ex.childs), make_function(std::forward<RHS>(rhs), placeholder<Is>()...));
  }

 // any other case e.g. f(x_+y_) = RHS etc .... which makes no sense : compiler will stop
 template<typename F, typename RHS, typename... T> 
  void operator<<(expr<tags::function, F, T...> && ex, RHS && rhs) = delete;
 template<typename F, typename RHS, typename... T> 
  void operator<<(expr<tags::function, F, T...> const & ex, RHS && rhs) = delete;

 /* ---------------------------------------------------------------------------------------------------
  * Auto assign for []
  *  --------------------------------------------------------------------------------------------------- */

 // by default it is deleted = not implemented : every class has to define it...
 template<typename T, typename F> void triqs_clef_auto_assign_subscript (T,F) = delete;

 // remove the ref_wrapper, terminal ...
 template<typename T, typename F> void triqs_clef_auto_assign_subscript (std::reference_wrapper<T> R     ,F && f) 
 { triqs_clef_auto_assign_subscript(R.get(),std::forward<F>(f));}
 template<typename T, typename F> void triqs_clef_auto_assign_subscript (expr<tags::terminal,T> const & t,F && f) 
 { triqs_clef_auto_assign_subscript(std::get<0>(t.childs),std::forward<F>(f));}

 // auto assign of an expr ? (for chain calls) : just reuse the same operator
 template<typename Tag, typename... Childs, typename RHS> 
  void triqs_clef_auto_assign_subscript (expr<Tag,Childs...> && ex, RHS const & rhs) { ex << rhs;}

 template<typename Tag, typename... Childs, typename RHS> 
  void triqs_clef_auto_assign_subscript (expr<Tag,Childs...> const & ex, RHS const & rhs) { ex << rhs;}

 // Same thing for the  [ ]
 template<typename F, typename RHS, int... Is> 
  void operator<<(expr<tags::subscript, F, placeholder<Is>...> const & ex, RHS && rhs) { 
   triqs_clef_auto_assign_subscript(std::get<0>(ex.childs), make_function(std::forward<RHS>(rhs), placeholder<Is>()...));
  }
  template<typename F, typename RHS, int... Is> 
  void operator<<(expr<tags::subscript, F, placeholder<Is>...> && ex, RHS && rhs) { 
   triqs_clef_auto_assign_subscript(std::get<0>(ex.childs), make_function(std::forward<RHS>(rhs), placeholder<Is>()...));
  }

  template<typename F, typename RHS, typename... T> 
  void operator<<(expr<tags::subscript, F, T...> && ex, RHS && rhs) = delete;

  template<typename F, typename RHS, typename... T> 
  void operator<<(expr<tags::subscript, F, T...> const & ex, RHS && rhs) = delete;

 /* --------------------------------------------------------------------------------------------------
  * Create a terminal node of an object. the from clone version force copying the object  
  * --------------------------------------------------------------------------------------------------- */

 // make a node with the ref, unless it is an rvalue (which is moved).
 template<typename T> expr<tags::terminal,typename expr_storage_t<T>::type >
  make_expr(T && x){ return {tags::terminal(), std::forward<T>(x)};}

 // make a node from a copy of the object
 template<typename T> expr<tags::terminal,typename remove_cv_ref<T>::type >
  make_expr_from_clone(T && x){ return {tags::terminal(), std::forward<T>(x)};}

 /* --------------------------------------------------------------------------------------------------
  * Create a call node of an object
  * The object can be kept as a : a ref, a copy, a view
  * --------------------------------------------------------------------------------------------------- */

 template<typename T> struct arity { static constexpr int value =-1;};

 namespace _result_of { 
  template< typename Obj, typename... Args > struct make_expr_call : 
   std::enable_if< is_any_lazy<Args...>::value, expr<tags::function,typename expr_storage_t<Obj>::type, typename expr_storage_t<Args>::type ...> > {
   static_assert (((arity<Obj>::value==-1) || (arity<Obj>::value == sizeof...(Args))), "Object called with a wrong number of arguments"); 
   };
 }
 template< typename Obj, typename... Args >
  typename _result_of::make_expr_call<Obj,Args...>::type 
  make_expr_call(Obj&& obj, Args &&... args) { return {tags::function(),std::forward<Obj>(obj), std::forward<Args>(args)...};}

 /* --------------------------------------------------------------------------------------------------
  * Create a [] call (subscript) node of an object
  * The object can be kept as a : a ref, a copy, a view
  * --------------------------------------------------------------------------------------------------- */

 namespace _result_of { 
  template< typename Obj, typename Arg> struct make_expr_subscript : 
   std::enable_if< is_any_lazy<Arg>::value, expr<tags::subscript,typename expr_storage_t<Obj>::type, typename expr_storage_t<Arg>::type> > {};
 }
 template< typename Obj, typename Arg>
  typename _result_of::make_expr_subscript<Obj,Arg>::type 
  make_expr_subscript(Obj&& obj, Arg && arg) { return {tags::subscript(),std::forward<Obj>(obj), std::forward<Arg>(arg)};}

 /* --------------------------------------------------------------------------------------------------
  *  function class : stores any expression polymorphically
  *  f(x_,y_ ) = an expression associates this expression dynamically to f, which 
  *  can then be used as a std::function of the same signature...
  * --------------------------------------------------------------------------------------------------- */
 template<typename F> class function;

 template<typename ReturnType, typename... T> class function<ReturnType(T...)> : tags::function_class  {
  typedef std::function<ReturnType(T...)> std_function_type;
  mutable std::shared_ptr <void>  _exp; // CLEAN THIS MUTABLE ?
  mutable std::shared_ptr < std_function_type > _fnt_ptr;
  public:
  function():_fnt_ptr{std::make_shared<std_function_type> ()}{}

  template<typename Expr, typename...  X>
   explicit function(Expr const & _e, X... x) : _exp(new Expr(_e)),_fnt_ptr(new std_function_type(make_function(_e, x...))){} 

  ReturnType operator()(T const &... t) const { return (*_fnt_ptr)(t...);}

#ifndef TRIQS_COMPILER_OBSOLETE_GCC
  template< typename... Args>
   auto operator()( Args&&... args ) const DECL_AND_RETURN(make_expr_call (*this,std::forward<Args>(args)...));
#else
  template< typename... Args>
  typename _result_of::make_expr_call<function const & ,Args...>::type
   operator()( Args&&... args ) const { return make_expr_call (*this,args...);}
#endif

  template<typename RHS> friend void triqs_clef_auto_assign (function const & x, RHS rhs) {
   * (x._fnt_ptr) =  std_function_type (rhs);
   x._exp = std::shared_ptr <void> (new typename std::remove_cv<decltype(rhs.ex)>::type (rhs.ex));
  }

 };
 template<typename F> struct force_copy_in_expr <function<F>> : std::true_type{};
  
 /* --------------------------------------------------------------------------------------------------
  *  The macro to make any function lazy
  *  TRIQS_CLEF_MAKE_FNT_LAZY (Arity,FunctionName ) : creates a new function in the triqs::lazy namespace 
  *  taking expressions (at least one argument has to be an expression) 
  *  The lookup happens by ADL, so IT MUST BE USED IN THE triqs::lazy namespace
  * --------------------------------------------------------------------------------------------------- */
#define TRIQS_CLEF_MAKE_FNT_LAZY(name)\
 struct name##_lazy_impl { \
  template<typename... A> auto operator()(A&&... a) const DECL_AND_RETURN (name(std::forward<A>(a)...));\
 };\
 template< typename... A> \
 auto name( A&& ... a) DECL_AND_RETURN(make_expr_call(name##_lazy_impl(),std::forward<A>(a)...));

#ifndef TRIQS_COMPILER_OBSOLETE_GCC

 #define TRIQS_CLEF_IMPLEMENT_LAZY_METHOD(TY,name)\
 struct __clef_lazy_method_impl_##name { \
  TY * _x;\
  template<typename... A> auto operator()(A&&... a) const DECL_AND_RETURN (_x->name(std::forward<A>(a)...));\
  friend std::ostream & operator<<(std::ostream & out, __clef_lazy_method_impl_##name  const & x) { return out<<BOOST_PP_STRINGIZE(TY)<<"."<<BOOST_PP_STRINGIZE(name);}\
 };\
 template< typename... A> \
 auto name( A&& ... a) DECL_AND_RETURN (make_expr_call(__clef_lazy_method_impl_##name{this},std::forward<A>(a)...));

#define TRIQS_CLEF_IMPLEMENT_LAZY_CALL(...)\
 template< typename... Args>\
 auto operator()(Args&&... args ) const & DECL_AND_RETURN(make_expr_call (*this,std::forward<Args>(args)...));\
 \
 template< typename... Args>\
 auto operator()(Args&&... args ) & DECL_AND_RETURN(make_expr_call (*this,std::forward<Args>(args)...));\
 \
 template< typename... Args>\
 auto operator()(Args&&... args ) && DECL_AND_RETURN(make_expr_call (std::move(*this),std::forward<Args>(args)...));\

#else

#define TRIQS_CLEF_IMPLEMENT_LAZY_METHOD(TY,name,RETURN_TYPE)\
 struct __clef_lazy_method_impl_##name { \
  TY * _x;\
  template<typename... A> RETURN_TYPE operator()(A&&... a) const {return ((*_x).name(std::forward<A>(a)...));}\
  friend std::ostream & operator<<(std::ostream & out, __clef_lazy_method_impl_##name  const & x) { return out<<BOOST_PP_STRINGIZE(TY)<<"."<<BOOST_PP_STRINGIZE(name);}\
 };\
 template< typename... A> \
 typename _result_of::make_expr_call<__clef_lazy_method_impl_##name, A...>::type\
 name( A&& ... a) {return make_expr_call(__clef_lazy_method_impl_##name{this},std::forward<A>(a)...);}

#define TRIQS_CLEF_IMPLEMENT_LAZY_CALL(TY)\
 template< typename... Args>\
 typename triqs::clef::_result_of::make_expr_call<TY const &,Args...>::type\
 operator()(Args&&... args ) const {return make_expr_call (*this,std::forward<Args>(args)...);}\
 \
 template< typename... Args>\
 typename triqs::clef::_result_of::make_expr_call<TY &,Args...>::type\
 operator()(Args&&... args ) {return make_expr_call (*this,std::forward<Args>(args)...);}

#endif

}} //  namespace triqs::clef
#endif 
