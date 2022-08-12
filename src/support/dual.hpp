//===-- duals/dual - Dual number class --------------------------*- C++ -*-===//
//
// Part of the cppduals project.
// https://tesch1.gitlab.io/cppduals
//
// (c)2019 Michael Tesch. tesch1@gmail.com
//
// See https://gitlab.com/tesch1/cppduals/blob/master/LICENSE.txt for
// license information.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef CPPDUALS_DUAL
#define CPPDUALS_DUAL

#ifndef PARSED_BY_DOXYGEN
#include <cmath>
#include <ctgmath>
#include <limits>
#include <type_traits>
#include <complex>
#include <random>
#include <iostream>
#endif

#if !defined(CPPDUALS_IGNORE_COMPILER_VERSION) && !defined(_WIN32)
#if __cplusplus < 201103L
#error CPPDUALS needs at least a C++11 compliant compiler
#endif
#endif

/// Configure whether system has POSIX extern int signgam;
#ifndef CPPDUALS_HAVE_SIGNGAM
#ifndef _WIN32
#define CPPDUALS_HAVE_SIGNGAM 1
#endif
#endif
namespace duals {

/**
\file       dual
\brief      Dual number class.

\mainpage   cppduals

\ref duals/dual is a single-header [Dual
number](https://en.wikipedia.org/wiki/Dual_number) C++ template
library that implements numbers of the type \f$(a + b \cdot
\epsilon)\f$ where \f$ \epsilon \ne 0\f$, and \f$\epsilon^2 = 0\f$.

`duals::dual<>` can be used for "automatic" differentiation, it can
also recursively nest with itself for higher orders of differentiation
(ie for the second derivative use `duals::dual<duals::dual<T>>`).  It
can also be used to differentiate parts of complex functions as
`std::complex<duals::dual<T>>`.  This file can be used stand-alone for
dual number support, but for Eigen vectorization support rather the
file \ref duals/dual_eigen should be included.

```
#include <duals/dual>

using namespace duals::literals;

template <class T> T   f(T x) { return pow(x,pow(x,x)); }
template <class T> T  df(T x) { return pow(x,-1. + x + pow(x,x)) * (1. + x*log(x) +
x*pow(log(x),2.)); } template <class T> T ddf(T x) { return (pow(x,pow(x,x)) * pow(pow(x,x - 1.) +
pow(x,x)*log(x)*(log(x) + 1.), 2.) + pow(x,pow(x,x)) * (pow(x,x - 1.) * log(x) + pow(x,x - 1.) *
(log(x) + 1.) + pow(x,x - 1.) * ((x - 1.)/x + log(x)) + pow(x,x) * log(x) * pow(log(x) + 1., 2.) ));
}

int main()
{
  std::cout << "  f(2.)            = " << f(2.)    << "\n";
  std::cout << " df(2.)            = " << df(2.)   << "\n";
  std::cout << "ddf(2.)            = " << ddf(2.)  << "\n";
  std::cout << "  f(2+1_e)         = " << f(2+1_e) << "\n";
  std::cout << "  f(2+1_e).dpart() = " << f(2+1_e).dpart() << "\n";
  duals::hyperduald x(2+1_e,1+0_e);
  std::cout << "  f((2+1_e) + (1+0_e)_e).dpart().dpart() = " << f(x).dpart().dpart() << "\n";
}
```
Produces (notice the derivative in the dual-part):
```
  f(2.)            = 16
 df(2.)            = 107.11
ddf(2.)            = 958.755
  f(2+1_e)         = (16+107.11_e)
  f(2+1_e).dpart() = 107.11
  f((2+1_e) + (1+0_e)_e).dpart().dpart() = 958.755
```

How this works can be seen by inspecting the infinite [Taylor
series](https://en.wikipedia.org/wiki/Taylor_series) expansion of a
function \f$ f(x) \f$ at \f$ a \f$ with \f$ x = a + b\epsilon \f$.
The series truncates itself due to the property \f$ \epsilon^2 = 0
\f$, leaving the function's derivative at \f$ a \f$ in the "dual part"
of the result (when \f$ b = 1 \f$):

\f[
\begin{split}
f(a + b \epsilon) &= f(a) + f'(a)(b \epsilon) + \frac{f''(a)}{2!}(b \epsilon)^2
+ \frac{f'''(a)}{3!}(b \epsilon)^3 + \ldots   \\        \
&= f(a) + f'(a)(b \epsilon)
\end{split}
\f]

The class is contained in a single c++11 header `#include
<duals/dual>`, for extended [Eigen](http://eigen.tuxfamily.org)
support, instead include the header \ref duals/dual_eigen "#include <duals/dual_eigen>".

Type X in the templates below can be any value which can be
assigned to value_type.

Type X also indicates a limitation to dual numbers of the same depth
but (possibly) different value_type as `duals::dual<T>`.  For example,
you can assign (or add/sub/mul/div) `duals::dual<float>` and
`duals::dual<double>`, but you can not assign
`duals::dual<duals::dual<float>>` to `duals::dual<float>`.

Here is a synopsis of the class:

```
namespace duals {

template<class T> class dual {

    typedef T value_type;

    dual(const & re = T(), const & du = T());
    dual(const dual &);
    template<class X> dual(const dual<X> &);

    T rpart() const;
    T dpart() const;

    void rpart(T);
    void dpart(T);

    dual<T> operator-() const;
    dual<T> operator+() const;

    dual<T> & operator= (const T &);
    dual<T> & operator+=(const T &);
    dual<T> & operator-=(const T &);
    dual<T> & operator*=(const T &);
    dual<T> & operator/=(const T &);

    dual<T> & operator=(const dual<T> &);
    template<class X> dual<T> & operator= (const dual<X> &);
    template<class X> dual<T> & operator+=(const dual<X> &);
    template<class X> dual<T> & operator-=(const dual<X> &);
    template<class X> dual<T> & operator*=(const dual<X> &);
    template<class X> dual<T> & operator/=(const dual<X> &);

    // The comparison operators are not strictly well-defined,
    // they are implemented as comparison of the real part.

    bool operator ==(const X &b) const;
    bool operator !=(const X &b) const;

    bool operator <(const X &b) const;
    bool operator >(const X &b) const;
    bool operator <=(const X &b) const;
    bool operator >=(const X &b) const;
    bool operator <(const dual<X> &b) const;
    bool operator >(const dual<X> &b) const;
    bool operator <=(const dual<X> &b) const;
    bool operator >=(const dual<X> &b) const;

};

// Non-member functions:

T rpart(dual<T>)        // Real part
T dpart(dual<T>)        // Dual part
dual<T> dconj(dual<T>)  // Dual-conjugate

dual<T> random(dual<T> low = {0,0}, dual<T> high = {1,0})

dual<T> exp(dual<T>)
dual<T> log(dual<T>)
dual<T> log10(dual<T>)
dual<T> pow(dual<T>, U)
dual<T> pow(U, dual<T>)
dual<T> pow(dual<T>, dual<T>)
dual<T> sqrt(dual<T>)
dual<T> cbrt(dual<T>)
dual<T> sin(dual<T>)
dual<T> cos(dual<T>)
dual<T> tan(dual<T>)
dual<T> asin(dual<T>)
dual<T> acos(dual<T>)
dual<T> atan(dual<T>)

// TODO:
dual<T> sinh(dual<T>)
dual<T> cosh(dual<T>)
dual<T> tanh(dual<T>)
dual<T> asinh(dual<T>)
dual<T> acosh(dual<T>)
dual<T> atanh(dual<T>)

// Non-differentiable operations on the real part.
T frexp(duals::dual<T> arg, int* exp );
duals::dual<T> ldexp(duals::dual<T> arg, int exp );
T trunc(duals::dual<T> d);
T floor(duals::dual<T> d);
T ceil(duals::dual<T> d);
T round(duals::dual<T> d);
int fpclassify(duals::dual<T> d);
bool isfinite(duals::dual<T> d);
bool isnormal(duals::dual<T> d);
bool isinf(duals::dual<T> d);
bool isnan(duals::dual<T> d);

// Stream IO
template<T> operator>>(basic_istream<charT, traits> &, dual<T> &);
template<T> operator<<(basic_ostream<charT, traits> &, const dual<T> &);

}
```

Some useful typedefs:

```
typedef dual<float> dualf;
typedef dual<double> duald;
typedef dual<long double> dualld;
typedef dual<dualf> hyperdualf;
typedef dual<duald> hyperduald;
typedef dual<dualld> hyperdualld;
typedef std::complex<dualf> cdualf;
typedef std::complex<duald> cduald;
typedef std::complex<dualld> cdualld;
```

There are also literals for dual parts defined in the inline
namespace duals::literals.  They are available with `using
namespace duals;` or `using namespace duals::literals`.

```
using namespace duals::literals;
dualf  x = 3 + 4_ef;
duald  y = 3 + 4_e;
dualld z = 3 + 4_el;
```

And in case you dislike iostreams, there are some formatters for the
[`{fmt}`](https://github.com/fmtlib/fmt) formatting library.  These
are disabled by default, but can be enabled by `#define
CPPDUALS_LIBFMT` for the `dual<>` formatter, and/or `#define
CPPDUALS_LIBFMT_COMPLEX` for the std::complex<> formatter. There are
three custom formatting flags that control how the dual numbers are
printed, these must come before the normal `{fmt}` formatting spec:
'$', '*', ','.  For me these are about 3x faster than iostreams.

```
#define CPPDUALS_LIBFMT
#define CPPDUALS_LIBFMT_COMPLEX
#inlcude <duals/dual>
using namespace duals::literals;
  ...
  string s = fmt::format("{:}",  1 + 2_e); // s = "(1.0+2.0_e)"
  string s = fmt::format("{:g}", 1 + 2_e); // s = "(1+2_e)"
  string s = fmt::format("{:*}", 1 + 2_e); // s = "(1.0+2.0*e)"
  string s = fmt::format("{:,}", 1 + 2_e); // s = "(1.0,2.0)"
  string s = fmt::format("{:*,g}", complexd(1,2_e)); // s = "((1)+(0,2)*i)"
```

The "archaic Greek epsilon" logo is from [Wikimedia
commons](https://commons.wikimedia.org/wiki/File:Greek_Epsilon_archaic.svg)

Some casual reading material:

- [ON QUATERNIONS, William Rowan Hamilton, Proceedings of the Royal Irish Academy, 3 (1847), pp.
1–16.](https://www.maths.tcd.ie/pub/HistMath/People/Hamilton/Quatern2/)
- [Basic Space-Time Transformations Expressed by Means of Two-Component Number
Systems](https://doi.org/10.12693%2Faphyspola.86.291)

 */

#ifndef PARSED_BY_DOXYGEN
template <class T> class dual;
#endif

/// Check if T is dual<>, match non-duals.
template <class T> struct is_dual : std::false_type {};

#ifndef PARSED_BY_DOXYGEN

/// Check if T is dual<>, match dual<>.
template <class T> struct is_dual<dual<T> > : std::true_type {};

#endif

/// Check if T is std::complex<>, match non- std::complex<>.
template <class T> struct is_complex : std::false_type {};

#ifndef PARSED_BY_DOXYGEN

/// Check if T is std::complex<>, match std::complex<>.
template <class T> struct is_complex<std::complex<T> > : std::true_type {};

#endif

/// Dual_traits helper class.
template <class T> struct dual_traits {
  /// Depth of T - for T=scalar this is 0. for dual_traits<double> it
  /// is 1.
  enum { depth = 0 }; // -Wenum-compare

  /// The real storage type.
  typedef T real_type;
};

#ifndef PARSED_BY_DOXYGEN

/// dual_traits for dual<> types
template <class T> struct dual_traits<dual<T> > {
  /// Depth to which this dual<> type is nested.  One (1) is a
  /// first-level dual, whereas non-duals have a depth of 0.
  enum { depth = dual_traits<T>::depth + 1 };

  /// The real storage type.
  typedef typename dual_traits<T>::real_type real_type;
};

template <class T> struct dual_traits<std::complex<dual<T> > > {
  /// complex<dual<T>> have the same 'depth' as their dual.
  enum { depth = dual_traits<T>::depth };

  /// The real storage type.
  typedef typename dual_traits<T>::real_type real_type;
};

namespace detail {

template <class T> struct Void { typedef void type; };
template <class T, class U = void> struct has_member_type : std::false_type {};
template <class T>
struct has_member_type<T, typename Void<typename T::type>::type> : std::true_type {
  struct wrap {
    typedef typename T::type type;
    typedef typename T::type ReturnType;
  };
};

} // namespace detail

/// Promote two types - default according to common_type
template <class T, class U, class V = void> struct promote : std::common_type<T, U> {};

/// Can types A and B be promoted to a common type?
template <class A, class B> using can_promote = detail::has_member_type<promote<A, B> >;

template <class T, class U>
struct promote<
    dual<T>, dual<U>,
    typename std::enable_if<(can_promote<T, U>::value &&
                             (int)dual_traits<T>::depth == (int)dual_traits<U>::depth)>::type> {
  typedef dual<typename promote<U, T>::type> type;
};
template <class T, class U>
struct promote<dual<T>, U,
               typename std::enable_if<(can_promote<T, U>::value &&
                                        (int)dual_traits<T>::depth >= (int)dual_traits<U>::depth &&
                                        !is_complex<U>::value)>::type> {
  typedef dual<typename promote<U, T>::type> type;
};
template <class T, class U>
struct promote<U, dual<T>,
               typename std::enable_if<(can_promote<T, U>::value &&
                                        (int)dual_traits<T>::depth >= (int)dual_traits<U>::depth &&
                                        !is_complex<U>::value)>::type> {
  typedef dual<typename promote<U, T>::type> type;
};
// /////////////////////////////////////////////////
template <class T, class U>
struct promote<std::complex<T>, std::complex<U>,
               typename std::enable_if<(can_promote<T, U>::value &&
                                        (is_dual<T>::value || is_dual<U>::value))>::type> {
  typedef std::complex<typename promote<T, U>::type> type;
};
template <class T, class U>
struct promote<
    std::complex<T>, U,
    typename std::enable_if<(can_promote<T, U>::value && (is_dual<T>::value || is_dual<U>::value) &&
                             !is_complex<U>::value)>::type> {
  typedef std::complex<typename promote<T, U>::type> type;
};
template <class T, class U>
struct promote<
    U, std::complex<T>,
    typename std::enable_if<(can_promote<T, U>::value && (is_dual<T>::value || is_dual<U>::value) &&
                             !is_complex<U>::value)>::type> {
  typedef std::complex<typename promote<T, U>::type> type;
};
// /////////////////////////////////////////////////

#endif // PARSED_BY_DOXYGEN

} // namespace duals

#ifndef PARSED_BY_DOXYGEN

#define NOMACRO // thwart macroification
#ifdef EIGEN_PI
#define MY_PI EIGEN_PI
#else
#define MY_PI M_PI
#endif

#endif // PARSED_BY_DOXYGEN

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_BEGIN_NAMESPACE_STD
#else
namespace std {
#endif

#ifdef CPPDUALS_ENABLE_STD_IS_ARITHMETIC

/// Duals are as arithmetic as their value_type is arithmetic.
template <class T> struct is_arithmetic<duals::dual<T> > : is_arithmetic<T> {};

#endif // CPPDUALS_ENABLE_IS_ARITHMETIC

/// Duals are compound types.
template <class T> struct is_compound<duals::dual<T> > : true_type {};

// Modification of std::numeric_limits<> per
// C++03 17.4.3.1/1, and C++11 18.3.2.3/1.
template <class T> struct numeric_limits<duals::dual<T> > : numeric_limits<T> {
  static constexpr bool is_specialized = true;
  static constexpr duals::dual<T> min NOMACRO() { return numeric_limits<T>::min NOMACRO(); }
  static constexpr duals::dual<T> lowest() { return numeric_limits<T>::lowest(); }
  static constexpr duals::dual<T> max NOMACRO() { return numeric_limits<T>::max NOMACRO(); }
  static constexpr duals::dual<T> epsilon() { return numeric_limits<T>::epsilon(); }
  static constexpr duals::dual<T> round_error() { return numeric_limits<T>::round_error(); }
  static constexpr duals::dual<T> infinity() { return numeric_limits<T>::infinity(); }
  static constexpr duals::dual<T> quiet_NaN() { return numeric_limits<T>::quiet_NaN(); }
  static constexpr duals::dual<T> signaling_NaN() { return numeric_limits<T>::signaling_NaN(); }
  static constexpr duals::dual<T> denorm_min() { return numeric_limits<T>::denorm_min(); }
};

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_END_NAMESPACE_STD
#else
} // namespace std
#endif

namespace duals {

#ifndef PARSED_BY_DOXYGEN

// T and X are wrapped in a dual<>
#define CPPDUALS_ONLY_SAME_DEPTH_AS_T(T, X)                                                        \
  typename std::enable_if<(int)duals::dual_traits<X>::depth == (int)duals::dual_traits<T>::depth,  \
                          int>::type = 0,                                                          \
                          typename std::enable_if<can_promote<T, X>::value, int>::type = 0

// Both T and U are wrapped in a dual<>
#define CPPDUALS_ENABLE_SAME_DEPTH_AND_COMMON_T(T, U)                                              \
  typename std::enable_if<(int)duals::dual_traits<T>::depth == (int)duals::dual_traits<U>::depth,  \
                          int>::type = 0,                                                          \
                          typename std::enable_if<can_promote<T, U>::value, int>::type = 0,        \
                          typename common_t = dual<typename duals::promote<T, U>::type>

// T is wrapped in a dual<>, U may or may not be.
#define CPPDUALS_ENABLE_LEQ_DEPTH_AND_COMMON_T(T, U)                                               \
  typename std::enable_if<((int)duals::dual_traits<T>::depth >=                                    \
                           (int)duals::dual_traits<U>::depth),                                     \
                          int>::type = 0,                                                          \
                          typename std::enable_if<can_promote<dual<T>, U>::value, int>::type = 0,  \
                          typename common_t = typename duals::promote<dual<T>, U>::type

// T is wrapped in complex<dual<>>
#define CPPDUALS_ENABLE_LEQ_DEPTH_AND_CX_COMMON_T(T, U)                                            \
  typename std::enable_if<                                                                         \
      ((int)duals::dual_traits<T>::depth >= (int)duals::dual_traits<U>::depth),                    \
      int>::type = 0,                                                                              \
      typename std::enable_if<can_promote<std::complex<dual<T> >, U>::value, int>::type = 0,       \
      typename common_t = typename duals::promote<std::complex<dual<T> >, U>::type

#define CPPDUALS_ENABLE_IF(...) typename std::enable_if<(__VA_ARGS__), int>::type = 0

#endif

/*! \page user Background
  TODO: Add text here...
*/

/// Abstract dual number class.  Can nest with other dual numbers and
/// complex numbers.
template <class T> class dual {
public:
  /// Class type of rpart() and dpart().  This type can be nested
  /// dual<> or std::complex<>.
  typedef T value_type;

private:
  /// The real part.
  value_type _real;
  /// The dual part.
  value_type _dual;

public:
  /// Construct dual from optional real and dual parts.
  constexpr dual(const value_type re = value_type(), const value_type du = value_type())
      : _real(re), _dual(du) {}

  /// Copy construct from a dual of equal depth.
  template <class X, CPPDUALS_ONLY_SAME_DEPTH_AS_T(T, X), CPPDUALS_ENABLE_IF(!is_complex<X>::value)>
  dual(const dual<X> &x) : _real(x.rpart()), _dual(x.dpart()) {}

  /// Cast to a complex<dual<>> with real part equal to *this.
  operator std::complex<dual<T> >() { return std::complex<dual<T> >(*this); }

  /// Explicit cast to an arithmetic type retains the rpart()
  template <class X, CPPDUALS_ENABLE_IF(std::is_arithmetic<X>::value && !is_dual<X>::value)>
  explicit operator X() const {
    return X(_real);
  }

  /// Get the real part.
  T rpart() const { return _real; }

  /// Get the dual part.
  T dpart() const { return _dual; }

  /// Set the real part.
  void rpart(value_type re) { _real = re; }

  /// Get the dual part.
  void dpart(value_type du) { _dual = du; }

  /// Unary negation
  dual<T> operator-() const { return dual<T>(-_real, -_dual); }

  /// Unary nothing
  dual<T> operator+() const { return *this; }

  /// Assignment of `value_type` assigns the real part and zeros the dual part.
  dual<T> &operator=(const T &x) {
    _real = x;
    _dual = value_type();
    return *this;
  }

  /// Add a relatively-scalar to this dual.
  dual<T> &operator+=(const T &x) {
    _real += x;
    return *this;
  }

  /// Subtract a relatively-scalar from this dual.
  dual<T> &operator-=(const T &x) {
    _real -= x;
    return *this;
  }

  /// Multiply a relatively-scalar with this dual.
  dual<T> &operator*=(const T &x) {
    _real *= x;
    _dual *= x;
    return *this;
  }

  /// Divide this dual by relatively-scalar.
  dual<T> &operator/=(const T &x) {
    _real /= x;
    _dual /= x;
    return *this;
  }

  // dua & operator=(const dual & x)  { _real =  x.rpart(); _dual =  x.dpart(); }
  template <class X, CPPDUALS_ONLY_SAME_DEPTH_AS_T(T, X)> dual<T> &operator=(const dual<X> &x) {
    _real = x.rpart();
    _dual = x.dpart();
    return *this;
  }

  /// Add a dual of the same depth to this dual.
  template <class X, CPPDUALS_ONLY_SAME_DEPTH_AS_T(T, X)> dual<T> &operator+=(const dual<X> &x) {
    _real += x.rpart();
    _dual += x.dpart();
    return *this;
  }

  /// Subtract a dual of the same depth from this dual.
  template <class X, CPPDUALS_ONLY_SAME_DEPTH_AS_T(T, X)> dual<T> &operator-=(const dual<X> &x) {
    _real -= x.rpart();
    _dual -= x.dpart();
    return *this;
  }

  /// Multiply this dual with a dual of same depth.
  template <class X, CPPDUALS_ONLY_SAME_DEPTH_AS_T(T, X)> dual<T> &operator*=(const dual<X> &x) {
    _dual = _real * x.dpart() + _dual * x.rpart();
    _real = _real * x.rpart();
    return *this;
  }

  /// Divide this dual by another dual of the same or lower depth.
  template <class X, CPPDUALS_ONLY_SAME_DEPTH_AS_T(T, X)> dual<T> &operator/=(const dual<X> &x) {
    _dual = (_dual * x.rpart() - _real * x.dpart()) / (x.rpart() * x.rpart());
    _real = _real / x.rpart();
    return *this;
  }

  // The following comparison operators are not strictly well-defined,
  // they are implemented as comparison of the real parts.
#if 0
  /// Compare this dual with another dual, comparing real parts for equality.
  template<class X, CPPDUALS_ONLY_SAME_DEPTH_AS_T(T,X)>
  bool operator ==(const dual<X> &b) const { return _real == b.rpart(); }

  /// Compare this dual with another dual, comparing real parts for inequality.
  template<class X, CPPDUALS_ONLY_SAME_DEPTH_AS_T(T,X)>
  bool operator !=(const dual<X> &b) const { return _real != b.rpart(); }

  /// Compare real part against real part of b
  template<class X, CPPDUALS_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,X)>
  bool operator <(const dual<X> &b) const { return _real < b.rpart(); }

  /// Compare real part against real part of b
  template<class X, CPPDUALS_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,X)>
  bool operator >(const dual<X> &b) const { return _real > b.rpart(); }

  /// Compare real part against real part of b
  template<class X, CPPDUALS_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,X)>
  bool operator <=(const dual<X> &b) const { return _real <= b.rpart(); }

  /// Compare real part against real part of b
  template<class X, CPPDUALS_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,X)>
  bool operator >=(const dual<X> &b) const { return _real >= b.rpart(); }
#endif
};

/// Get the dual's real part.
template <class T> T rpart(const dual<T> &x) { return x.rpart(); }

/// Get the dual's dual part.
template <class T> T dpart(const dual<T> &x) { return x.dpart(); }

/// R-part of complex<dual<>> is non-dual complex<> (not to be confused with real())
template <class T> std::complex<T> rpart(const std::complex<dual<T> > &x) {
  return std::complex<T>(x.real().rpart(), x.imag().rpart());
}

/// Dual part of complex<dual<>> is complex<>
template <class T> std::complex<T> dpart(const std::complex<dual<T> > &x) {
  return std::complex<T>(x.real().dpart(), x.imag().dpart());
}

/// Get a non-dual's real part.
template <class T,
          CPPDUALS_ENABLE_IF((std::is_arithmetic<T>::value && !std::is_compound<T>::value) ||
                             is_complex<T>::value)>
T rpart(const T &x) {
  return x;
}

/// Get a non-dual's dual part := zero.
template <class T,
          CPPDUALS_ENABLE_IF((std::is_arithmetic<T>::value && !std::is_compound<T>::value) ||
                             is_complex<T>::value)>
T dpart(const T &) {
  return T(0);
}

#ifndef PARSED_BY_DOXYGEN

/// Dual +-*/ ops with another entity
#define CPPDUALS_BINARY_OP(op)                                                                     \
  template <class T, class U, CPPDUALS_ENABLE_SAME_DEPTH_AND_COMMON_T(T, U)>                       \
  common_t operator op(const dual<T> &z, const dual<U> &w) {                                       \
    common_t x(z);                                                                                 \
    return x op## = w;                                                                             \
  }                                                                                                \
  template <class T, class U, CPPDUALS_ENABLE_LEQ_DEPTH_AND_COMMON_T(T, U),                        \
            CPPDUALS_ENABLE_IF(!std::is_same<U, std::complex<dual<T> > >::value)>                  \
  common_t operator op(const dual<T> &z, const U &w) {                                             \
    common_t x(z);                                                                                 \
    return x op## = w;                                                                             \
  }                                                                                                \
  template <class T, class U, CPPDUALS_ENABLE_LEQ_DEPTH_AND_COMMON_T(T, U),                        \
            CPPDUALS_ENABLE_IF(!std::is_same<U, std::complex<dual<T> > >::value)>                  \
  common_t operator op(const U &z, const dual<T> &w) {                                             \
    common_t x(z);                                                                                 \
    return x op## = w;                                                                             \
  }                                                                                                \
  template <class T, class U, CPPDUALS_ENABLE_LEQ_DEPTH_AND_CX_COMMON_T(T, U),                     \
            CPPDUALS_ENABLE_IF(!std::is_same<U, std::complex<dual<T> > >::value)>                  \
  common_t operator op(const std::complex<dual<T> > &z, const U &w) {                              \
    common_t x(z);                                                                                 \
    return x op## = w;                                                                             \
  }                                                                                                \
  template <class T, class U, CPPDUALS_ENABLE_LEQ_DEPTH_AND_CX_COMMON_T(T, U),                     \
            CPPDUALS_ENABLE_IF(!std::is_same<U, std::complex<dual<T> > >::value)>                  \
  common_t operator op(const U &z, const std::complex<dual<T> > &w) {                              \
    common_t x(z);                                                                                 \
    return x op## = w;                                                                             \
  }

CPPDUALS_BINARY_OP(+)
CPPDUALS_BINARY_OP(-)
CPPDUALS_BINARY_OP(*)
CPPDUALS_BINARY_OP(/)

/// Dual compared to a non-complex lower rank thing
#define CPPDUALS_LHS_COMPARISON(op)                                                                \
  template <class T, class U, CPPDUALS_ENABLE_SAME_DEPTH_AND_COMMON_T(T, U)>                       \
  bool operator op(const dual<T> &a, const dual<U> &b) {                                           \
    return a.rpart() op b.rpart();                                                                 \
  }                                                                                                \
  template <class T, class U, CPPDUALS_ENABLE_LEQ_DEPTH_AND_COMMON_T(T, U),                        \
            CPPDUALS_ENABLE_IF(!is_complex<U>::value)>                                             \
  bool operator op(const U &a, const dual<T> &b) {                                                 \
    return a op b.rpart();                                                                         \
  }                                                                                                \
  template <class T, class U, CPPDUALS_ENABLE_LEQ_DEPTH_AND_COMMON_T(T, U),                        \
            CPPDUALS_ENABLE_IF(!is_complex<U>::value)>                                             \
  bool operator op(const dual<T> &a, const U &b) {                                                 \
    return a.rpart() op b;                                                                         \
  }

CPPDUALS_LHS_COMPARISON(<)
CPPDUALS_LHS_COMPARISON(>)
CPPDUALS_LHS_COMPARISON(<=)
CPPDUALS_LHS_COMPARISON(>=)
CPPDUALS_LHS_COMPARISON(==)
CPPDUALS_LHS_COMPARISON(!=)

#endif // PARSED_BY_DOXYGEN

} // namespace duals

// some useful functions of solely the real part
#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_BEGIN_NAMESPACE_STD
#else
namespace std {
#endif

#ifndef PARSED_BY_DOXYGEN

#define make_math(T)                                                                               \
  inline T(frexp)(const duals::dual<T> &arg, int *exp) { return (frexp)(arg.rpart(), exp); }       \
  inline duals::dual<T>(ldexp)(const duals::dual<T> &arg, int exp) {                               \
    return arg * std::pow((T)2, exp);                                                              \
  }                                                                                                \
  inline T(trunc)(const duals::dual<T> &d) { return (trunc)(d.rpart()); }                          \
  inline T(floor)(const duals::dual<T> &d) { return (floor)(d.rpart()); }                          \
  inline T(ceil)(const duals::dual<T> &d) { return (ceil)(d.rpart()); }                            \
  inline T(round)(const duals::dual<T> &d) { return (round)(d.rpart()); }                          \
  inline int(fpclassify)(const duals::dual<T> &d) { return (fpclassify)(d.rpart()); }              \
  inline bool(isfinite)(const duals::dual<T> &d) { return (isfinite)(d.rpart()); }                 \
  inline bool(isnormal)(const duals::dual<T> &d) { return (isnormal)(d.rpart()); }                 \
  inline bool(isinf)(const duals::dual<T> &d) { return (isinf)(d.rpart()); }                       \
  inline bool(isnan)(const duals::dual<T> &d) { return (isnan)(d.rpart()); }

make_math(float) make_math(double) make_math(long double)

#undef make_math

#endif // PARSED_BY_DOXYGEN

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
    _LIBCPP_END_NAMESPACE_STD
#else
} // namespace std
#endif
    namespace duals {

  namespace randos {

  // Random real value between a and b.
  template <class T, typename dist = std::uniform_real_distribution<T>,
            CPPDUALS_ENABLE_IF(!is_complex<T>::value && !is_dual<T>::value)>
  T random(T a = T(0), T b = T(1)) {
    static std::default_random_engine generator;
    dist distribution(a, b);
    return distribution(generator);
  }

  template <class T, typename dist = std::uniform_real_distribution<T>,
            CPPDUALS_ENABLE_IF(!is_complex<T>::value && !is_dual<T>::value)>
  T random2(T a = T(0), T b = T(1)) {
    return random<T, dist>(a, b);
  }

  // Helper class for testing - also random value in dual part.
  template <class DT, CPPDUALS_ENABLE_IF(is_dual<DT>::value)>
  DT random2(DT a = DT(0, 0), DT b = DT(1, 1)) {
    using randos::random;
    return DT(a.rpart() + random2<typename DT::value_type>() * (b.rpart() - a.rpart()),
              a.dpart() + random2<typename DT::value_type>() * (b.dpart() - a.dpart()));
  }

  // Helper class for testing - also random value in dual part of the complex.
  template <class CT, CPPDUALS_ENABLE_IF(is_complex<CT>::value)>
  CT random2(CT a = CT(0, 0), CT b = CT(1, 1)) {
    return CT(a.real() + random2<typename CT::value_type>() * (b.real() - a.real()),
              a.imag() + random2<typename CT::value_type>() * (b.imag() - a.imag()));
  }

  } // namespace randos

  /// Random real and dual parts, used by Eigen's Random(), by default
  // the returned value has zero dual part and is in the range [0+0_e,
  // 1+0_e].
  template <class DT, typename dist = std::uniform_real_distribution<typename DT::value_type>,
            CPPDUALS_ENABLE_IF(is_dual<DT>::value)>
  DT random(DT a = DT(0, 0), DT b = DT(1, 0)) {
    using randos::random;
    return DT(random<typename DT::value_type, dist>(a.rpart(), b.rpart()),
              random<typename DT::value_type, dist>(a.dpart(), b.dpart()));
  }

  /// Complex Conjugate of a dual is just the dual.
  template <class T> dual<T> conj(const dual<T> &x) { return x; }

  /// Conjugate a thing that's not dual and not complex -- it has no
  /// complex value so just return it.  This is different from
  /// std::conj() which promotes the T to a std::complex<T>.
  template <class T, CPPDUALS_ENABLE_IF(!is_dual<T>::value && !is_complex<T>::value &&
                                        std::is_arithmetic<T>::value)>
  T conj(const T &x) {
    return x;
  }

  /// Dual Conjugate, such that x*dconj(x) = rpart(x)^2.  Different
  /// function name from complex conjugate conj().
  template <class T> dual<T> dconj(const dual<T> &x) { return dual<T>(x.rpart(), -x.dpart()); }

  /// Dual Conjugate of a complex, such that x*dconj(x) = rpart(x)^2.
  /// Different function name from complex conjugate conj().
  template <class T> std::complex<T> dconj(const std::complex<T> &x) {
    return std::complex<T>(dconj(x.real()), dconj(x.imag()));
  }

  /// Conjugate a thing that's not dual and not complex.
  template <class T, CPPDUALS_ENABLE_IF(!is_dual<T>::value && !is_complex<T>::value &&
                                        std::is_arithmetic<T>::value)>
  T dconj(const T &x) {
    return x;
  }

  /// Exponential e^x
  template <class T> dual<T> exp(const dual<T> &x) {
    using std::exp;
    T v = exp(x.rpart());
    return dual<T>(v, v * x.dpart());
  }

  /// Natural log ln(x)
  template <class T> dual<T> log(const dual<T> &x) {
    using std::log;
    T v = log(x.rpart());
    if (x.dpart() == T(0))
      return v;
    else
      return dual<T>(v, x.dpart() / x.rpart());
  }

  template <class T> dual<T> log10(const dual<T> &x) {
    using std::log;
    return log(x) / log(static_cast<T>(10));
  }

  template <class T> dual<T> log2(const dual<T> &x) {
    using std::log;
    return log(x) / log(static_cast<T>(2));
  }

  template <class T, class U, CPPDUALS_ENABLE_LEQ_DEPTH_AND_COMMON_T(T, U)>
  common_t pow(const dual<T> &x, const U &y) {
    using std::pow;
    typedef typename common_t::value_type V;

    return common_t(pow(x.rpart(), y),
                    x.dpart() == V(0) ? V(0) : x.dpart() * y * pow(x.rpart(), y - U(1)));
  }

  template <class T, class U, CPPDUALS_ENABLE_LEQ_DEPTH_AND_COMMON_T(T, U)>
  common_t pow(const U &x, const dual<T> &y) {
    return pow(common_t(x), y);
  }

  template <class T, class U, CPPDUALS_ENABLE_SAME_DEPTH_AND_COMMON_T(T, U)>
  common_t pow(const dual<T> &f, const dual<U> &g) {
    using std::log;
    using std::pow;
#if 1
    using std::floor;
    typedef typename common_t::value_type V;
    common_t result;

    if (f.rpart() == T(0) && g.rpart() >= U(1)) {
      if (g.rpart() > U(1)) { result = common_t(0); }
      else { result = f; }
    }
    else {
      if (f.rpart() < T(0) && g.rpart() == floor(g.rpart())) {
        V const tmp = g.rpart() * pow(f.rpart(), g.rpart() - U(1.0));
        result = common_t(pow(f.rpart(), g.rpart()), f.dpart() == T(0) ? T(0) : f.dpart() * tmp);
        if (g.dpart() != U(0.0)) {
          // Return a NaN when g.dpart() != 0.
          result.dpart(std::numeric_limits<T>::quiet_NaN());
        }
      }
      else {
        // Handle the remaining cases. For cases 4,5,6,9 we allow the log()
        // function to generate -HUGE_VAL or NaN, since those cases result in a
        // nonfinite derivative.
        V const tmp1 = pow(f.rpart(), g.rpart());
        V const tmp2 = g.rpart() * pow(f.rpart(), g.rpart() - T(1.0));
        V const tmp3 = tmp1 * log(f.rpart());
        result = common_t(tmp1, f.dpart() == T(0) && g.dpart() == U(0)
                                    ? V(0)
                                    : tmp2 * f.dpart() + tmp3 * g.dpart());
      }
    }

    return result;
#else
  T v = pow(f.rpart(), g.rpart());
  return common_t(v, pow(f.rpart(), g.rpart() - T(1)) *
                         (g.rpart() * f.dpart() + f.rpart() * log(f.rpart()) * g.dpart()));
#endif
  }

  namespace utils {
  template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }
  } // namespace utils

  template <class T> dual<T> abs(const dual<T> &x) {
    using std::abs;
    return dual<T>(abs(x.rpart()), x.dpart() * utils::sgn(x.rpart()));
  }

  template <class T> dual<T> fabs(const dual<T> &x) {
    using std::fabs;
    return dual<T>(fabs(x.rpart()), x.dpart() * utils::sgn(x.rpart()));
  }

#if 0
// TODO?
template<class T> dual<T> abs2(const dual<T> & x) {
  using std::abs;
  return dual<T>(x.rpart() * x.rpart(),
                 xxx x.dpart() * utils::sgn(x.rpart()));
}
#endif

  template <class T> duals::dual<T> copysign(const duals::dual<T> &x, const duals::dual<T> &y) {
    using std::copysign;
    T r = copysign(x.rpart(), y.rpart());
    return duals::dual<T>(r, r == x.rpart() ? x.dpart() : -x.dpart());
  }

  template <class T> duals::dual<T> hypot(const duals::dual<T> &x, const duals::dual<T> &y) {
    return sqrt(x * x + y * y);
  }

  template <class T> duals::dual<T> scalbn(const duals::dual<T> &arg, int ex) {
    return arg * std::pow((T)2, ex);
  }

  template <class T> duals::dual<T>(fmax)(const duals::dual<T> &x, const duals::dual<T> &y) {
    return x.rpart() > y.rpart() ? x : y;
  }

  template <class T> duals::dual<T>(fmin)(const duals::dual<T> &x, const duals::dual<T> &y) {
    return x.rpart() <= y.rpart() ? x : y;
  }

  template <class T> duals::dual<T> logb(const duals::dual<T> &x) { return duals::log2(x); }

  template <class T> int(fpclassify)(const duals::dual<T> &d) {
    using std::fpclassify;
    return (fpclassify)(d.rpart());
  }
  template <class T> bool(isfinite)(const duals::dual<T> &d) {
    using std::isfinite;
    return (isfinite)(d.rpart());
  }
  template <class T> bool(isnormal)(const duals::dual<T> &d) {
    using std::isnormal;
    return (isnormal)(d.rpart());
  }
  template <class T> bool(isinf)(const duals::dual<T> &d) {
    using std::isinf;
    return (isinf)(d.rpart());
  }
  template <class T> bool(isnan)(const duals::dual<T> &d) {
    using std::isnan;
    return (isnan)(d.rpart());
  }
  template <class T> bool(signbit)(const duals::dual<T> &d) {
    using std::signbit;
    return (signbit)(d.rpart());
  }

  template <class T> dual<T> sqrt(const dual<T> &x) {
    using std::sqrt;
    T v = sqrt(x.rpart());
    if (x.dpart() == T(0))
      return v;
    else
      return dual<T>(v, x.dpart() / (T(2) * v));
  }

  template <class T> dual<T> cbrt(const dual<T> &x) {
    using std::cbrt;
    T v = cbrt(x.rpart());
    if (x.dpart() == T(0))
      return v;
    else
      return dual<T>(v, x.dpart() / (T(3) * v * v));
  }

  template <class T> dual<T> sin(const dual<T> &x) {
    using std::cos;
    using std::sin;
    return dual<T>(sin(x.rpart()), x.dpart() * cos(x.rpart()));
  }

  template <class T> dual<T> cos(const dual<T> &x) {
    using std::cos;
    using std::sin;
    return dual<T>(cos(x.rpart()), -sin(x.rpart()) * x.dpart());
  }

  template <class T> dual<T> tan(const dual<T> &x) {
    using std::tan;
    T v = tan(x.rpart());
    return dual<T>(v, x.dpart() * (v * v + 1));
  }

  template <class T> dual<T> asin(const dual<T> &x) {
    using std::asin;
    using std::sqrt;
    T v = asin(x.rpart());
    if (x.dpart() == T(0))
      return v;
    else
      return dual<T>(v, x.dpart() / sqrt(1 - x.rpart() * x.rpart()));
  }

  template <class T> dual<T> acos(const dual<T> &x) {
    using std::acos;
    using std::sqrt;
    T v = acos(x.rpart());
    if (x.dpart() == T(0))
      return v;
    else
      return dual<T>(v, -x.dpart() / sqrt(1 - x.rpart() * x.rpart()));
  }

  template <class T> dual<T> atan(const dual<T> &x) {
    using std::atan;
    T v = atan(x.rpart());
    if (x.dpart() == T(0))
      return v;
    else
      return dual<T>(v, x.dpart() / (1 + x.rpart() * x.rpart()));
  }

  template <class T> dual<T> atan2(const dual<T> &y, const dual<T> &x) {
    using std::atan2;
    T v = atan2(y.rpart(), x.rpart());
    return dual<T>(v, (x.rpart() * y.dpart() - y.rpart() * x.dpart()) /
                          (y.rpart() * y.rpart() + x.rpart() * x.rpart()));
  }

  // TODO
  template <class T, class U, CPPDUALS_ENABLE_LEQ_DEPTH_AND_COMMON_T(T, U)>
  common_t atan2(const dual<T> &y, const U &x) {
    using std::atan2;
    T v = atan2(y.rpart(), x);
    return dual<T>(v, (x * y.dpart()) / (y.rpart() * y.rpart() + x * x));
  }

  // TODO
  template <class T, class U, CPPDUALS_ENABLE_LEQ_DEPTH_AND_COMMON_T(T, U)>
  common_t atan2(const U &y, const dual<T> &x) {
    using std::atan2;
    T v = atan2(y, x.rpart());
    return dual<T>(v, (-y * x.dpart()) / (y * y + x.rpart() * x.rpart()));
  }

  // TODO
  template <class T> dual<T> sinh(const dual<T> &x);
  template <class T> dual<T> cosh(const dual<T> &x);
  template <class T> dual<T> tanh(const dual<T> &x);
  template <class T> dual<T> asinh(const dual<T> &x);
  template <class T> dual<T> acosh(const dual<T> &x);
  template <class T> dual<T> atanh(const dual<T> &x);
  template <class T> dual<T> log1p(const dual<T> &x);
  template <class T> dual<T> expm1(const dual<T> &x);

  ////////////////////////////////////////////////////////////////
  // added for meep forward-diff operations
  ////////////////////////////////////////////////////////////////

  /// The tanh function.  Make sure to `#include <math.h>` before
  /// `#include <duals/dual>` to use this function.
  template <class T> dual<T> tanh(const dual<T> &x) {
    using std::tanh;
    T sech = 1.0 / std::cosh(x.rpart());
    return dual<T>(tanh(x.rpart()), x.dpart() * sech * sech);
  }

  ////////////////////////////////////////////////////////////////
  //
  ////////////////////////////////////////////////////////////////

  /// The error function.  Make sure to `#include <math.h>` before
  /// `#include <duals/dual>` to use this function.
  template <class T> dual<T> erf(const dual<T> &x) {
    using std::erf;
    using std::exp;
    using std::pow;
    using std::sqrt;
    return dual<T>(erf(x.rpart()), x.dpart() * T(2) / sqrt(T(MY_PI)) * exp(-pow(x.rpart(), T(2))));
  }

  /// Error function complement (1 - erf()).
  template <class T> dual<T> erfc(const dual<T> &x) {
    using std::erfc;
    using std::exp;
    using std::pow;
    using std::sqrt;
    return dual<T>(erfc(x.rpart()),
                   x.dpart() * -T(2) / sqrt(T(MY_PI)) * exp(-pow(x.rpart(), T(2))));
  }

  /// Gamma function.  Approximation of the dual part.
  // TODO specialize for integers
  template <class T> dual<T> tgamma(const dual<T> &x) {
    using std::tgamma;
    T v = tgamma(x.rpart());
    if (x.dpart() == T(0))
      return v;
    else {
      int errno_saved = errno;
      T h(T(1) / (1ull << (std::numeric_limits<T>::digits / 3)));
      T w((tgamma(x.rpart() + h) - tgamma(x.rpart() - h)) / (2 * h));
      errno = errno_saved;
      return dual<T>(v, x.dpart() * w);
    }
  }

  /// Log of absolute value of gamma function.  Approximation of the dual part.
  template <class T> dual<T> lgamma(const dual<T> &x) {
    using std::lgamma;
    T v = lgamma(x.rpart());
    if (x.dpart() == T(0))
      return v;
    else {
#if CPPDUALS_HAVE_SIGNGAM
      int signgam_saved = signgam;
#endif
      int errno_saved = errno;
      T h(T(1) / (1ull << (std::numeric_limits<T>::digits / 3)));
      T w((lgamma(x.rpart() + h) - lgamma(x.rpart() - h)) / (2 * h));
#if CPPDUALS_HAVE_SIGNGAM
      signgam = signgam_saved;
#endif
      errno = errno_saved;
      return dual<T>(v, x.dpart() * w);
    }
  }

  /// Putto operator
  template <class T, class _CharT, class _Traits>
  std::basic_ostream<_CharT, _Traits> &operator<<(std::basic_ostream<_CharT, _Traits> &os,
                                                  const dual<T> &x) {
    std::basic_ostringstream<_CharT, _Traits> s;
    s.flags(os.flags());
    s.imbue(os.getloc());
    s.precision(os.precision());
    s << '(' << x.rpart() << (x.dpart() < 0 ? "" : "+") << x.dpart() << "_e"
      << (std::is_same<typename std::decay<T>::type, float>::value         ? "f"
          : std::is_same<typename std::decay<T>::type, double>::value      ? ""
          : std::is_same<typename std::decay<T>::type, long double>::value ? "l"
                                                                           : "")
      << ")";
    return os << s.str();
  }

  /// Stream reader
  template <class T, class CharT, class Traits>
  std::basic_istream<CharT, Traits> &operator>>(std::basic_istream<CharT, Traits> &is, dual<T> &x) {
    if (is.good()) {
      ws(is);
      if (is.peek() == CharT('(')) {
        is.get();
        T r;
        is >> r;
        if (!is.fail()) {
          CharT c = is.peek();
          if (c != CharT('_')) {
            ws(is);
            c = is.peek();
          }
          if (c == CharT('+') || c == CharT('-') || c == CharT('_')) {
            if (c == CharT('+')) is.get();
            T d;
            if (c == CharT('_')) {
              d = r;
              r = 0;
            }
            else
              is >> d;
            if (!is.fail()) {
              ws(is);
              c = is.peek();
              if (c == CharT('_')) {
                is.get();
                c = is.peek();
                if (c == CharT('e')) {
                  is.get();
                  c = is.peek();
                  if ((c == 'f' && !std::is_same<typename std::decay<T>::type, float>::value) ||
                      (c == 'l' && !std::is_same<typename std::decay<T>::type, long double>::value))
                    is.setstate(std::ios_base::failbit);
                  else {
                    if (c == 'f' || c == 'l') is.get();
                    ws(is);
                    c = is.peek();
                    if (c == ')') {
                      is.get();
                      x = dual<T>(r, d);
                    }
                    else
                      is.setstate(std::ios_base::failbit);
                  }
                }
                else
                  is.setstate(std::ios_base::failbit);
              }
              else
                is.setstate(std::ios_base::failbit);
            }
            else
              is.setstate(std::ios_base::failbit);
          }
          else if (c == CharT(')')) {
            is.get();
            x = dual<T>(r, T(0));
          }
          else
            is.setstate(std::ios_base::failbit);
        }
        else
          is.setstate(std::ios_base::failbit);
      }
      else {
        T r;
        is >> r;
        if (!is.fail())
          x = dual<T>(r, T(0));
        else
          is.setstate(std::ios_base::failbit);
      }
    }
    else
      is.setstate(std::ios_base::failbit);
    return is;
  }

#if __cpp_user_defined_literals >= 200809
  /// Dual number literals in namespace duals::literals
  inline namespace literals {
  using duals::dual;
  constexpr dual<float> operator"" _ef(long double du) { return {0.0f, static_cast<float>(du)}; }
  constexpr dual<float> operator"" _ef(unsigned long long du) {
    return {0.0f, static_cast<float>(du)};
  }
  constexpr dual<double> operator"" _e(long double du) { return {0.0, static_cast<double>(du)}; }
  constexpr dual<double> operator"" _e(unsigned long long du) {
    return {0.0, static_cast<double>(du)};
  }
  constexpr dual<long double> operator"" _el(long double du) { return {0.0l, du}; }
  constexpr dual<long double> operator"" _el(unsigned long long du) {
    return {0.0l, static_cast<long double>(du)};
  }
  } // namespace literals
#endif

  typedef dual<float> dualf;
  typedef dual<double> duald;
  typedef dual<long double> dualld;
  typedef dual<dualf> hyperdualf;
  typedef dual<duald> hyperduald;
  typedef dual<dualld> hyperdualld;
  typedef std::complex<dualf> cdualf;
  typedef std::complex<duald> cduald;
  typedef std::complex<dualld> cdualld;

} // namespace duals

#ifdef CPPDUALS_LIBFMT
#include <fmt/format.h>
/// duals::dual<> Formatter for libfmt https://github.com/fmtlib/fmt
///
/// Formats a dual number (r,d) as (r+d_e), offering the same
/// formatting options as the underlying type - with the addition of
/// three optional format options, only one of which may appear
/// directly after the ':' in the format spec: '$', '*', and ','.  The
/// '*' flag changes the separating _ to a *, producing (r+d*e), where
/// r and d are the formatted value_type values.  The ',' flag simply
/// prints the real and dual parts separated by a comma.  As a
/// concrete exmple, this formatter can produce either (3+5.4_e) or
/// (3+5.4*e) or (3,5.4) for a dual<double> using the specs {:g},
/// {:*g}, or {:,g}, respectively.  When the '*' is NOT specified, the
/// output should be compatible with the input operator>> and the dual
/// literals below. (this implementation is a bit hacky - glad for
/// cleanups).
template <typename T, typename Char>
struct fmt::formatter<duals::dual<T>, Char> : public fmt::formatter<T, Char> {
  typedef fmt::formatter<T, Char> base;
  enum style { expr, star, pair } style_ = expr;
  fmt::detail::dynamic_format_specs<Char> specs_;
  FMT_CONSTEXPR auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
    using handler_type = fmt::detail::dynamic_specs_handler<format_parse_context>;
    auto type = fmt::detail::type_constant<T, Char>::value;
    fmt::detail::specs_checker<handler_type> handler(handler_type(specs_, ctx), type);
    auto it = ctx.begin();
    if (it != ctx.end()) switch (*it) {
        case '$':
          style_ = style::expr;
          ctx.advance_to(++it);
          break;
        case '*':
          style_ = style::star;
          ctx.advance_to(++it);
          break;
        case ',':
          style_ = style::pair;
          ctx.advance_to(++it);
          break;
        default: break;
      }
    parse_format_specs(ctx.begin(), ctx.end(), handler);
    return base::parse(ctx);
  }
  template <typename FormatCtx>
  auto format(const duals::dual<T> &x, FormatCtx &ctx) -> decltype(ctx.out()) {
    format_to(ctx.out(), "(");
    if (style_ == style::pair) {
      base::format(x.rpart(), ctx);
      format_to(ctx.out(), ",");
      base::format(x.dpart(), ctx);
      return format_to(ctx.out(), ")");
    }
    if (x.rpart() || !x.dpart()) base::format(x.rpart(), ctx);
    if (x.dpart()) {
      if (x.rpart() && x.dpart() >= 0 && specs_.sign != sign::plus) format_to(ctx.out(), "+");
      base::format(x.dpart(), ctx);
      if (style_ == style::star)
        format_to(ctx.out(), "*e");
      else
        format_to(ctx.out(), "_e");
      if (std::is_same<typename std::decay<T>::type, float>::value) format_to(ctx.out(), "f");
      if (std::is_same<typename std::decay<T>::type, long double>::value) format_to(ctx.out(), "l");
    }
    return format_to(ctx.out(), ")");
  }
};
#endif

#ifdef CPPDUALS_LIBFMT_COMPLEX
#ifndef CPPDUALS_LIBFMT
#include <fmt/format.h>
#endif
/// std::complex<> Formatter for libfmt https://github.com/fmtlib/fmt
///
/// libfmt does not provide a formatter for std::complex<>, although
/// one is proposed for c++20.  Anyway, at the expense of a k or two,
/// you can define CPPDUALS_LIBFMT_COMPLEX and get this one.
///
/// The standard iostreams formatting of complex numbers is (a,b),
/// where a and b are the real and imaginary parts.  This formats a
/// complex number (a+bi) as (a+bi), offering the same formatting
/// options as the underlying type - with the addition of three
/// optional format options, only one of which may appear directly
/// after the ':' in the format spec (before any fill or align): '$'
/// (the default if no flag is specified), '*', and ','.  The '*' flag
/// adds a * before the 'i', producing (a+b*i), where a and b are the
/// formatted value_type values.  The ',' flag simply prints the real
/// and complex parts separated by a comma (same as iostreams' format).
/// As a concrete exmple, this formatter can produce either (3+5.4i)
/// or (3+5.4*i) or (3,5.4) for a complex<double> using the specs {:g}
/// | {:$g}, {:*g}, or {:,g}, respectively.  (this implementation is a
/// bit hacky - glad for cleanups).
///
template <typename T, typename Char>
struct fmt::formatter<std::complex<T>, Char> : public fmt::formatter<T, Char> {
  typedef fmt::formatter<T, Char> base;
  enum style { expr, star, pair } style_ = expr;
  fmt::detail::dynamic_format_specs<Char> specs_;
  FMT_CONSTEXPR auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
    using handler_type = fmt::detail::dynamic_specs_handler<format_parse_context>;
    auto type = fmt::detail::type_constant<T, Char>::value;
    fmt::detail::specs_checker<handler_type> handler(handler_type(specs_, ctx), type);
    auto it = ctx.begin();
    if (it != ctx.end()) {
      switch (*it) {
        case '$':
          style_ = style::expr;
          ctx.advance_to(++it);
          break;
        case '*':
          style_ = style::star;
          ctx.advance_to(++it);
          break;
        case ',':
          style_ = style::pair;
          ctx.advance_to(++it);
          break;
        default: break;
      }
    }
    parse_format_specs(ctx.begin(), ctx.end(), handler);
    // todo: fixup alignment
    return base::parse(ctx);
  }
  template <typename FormatCtx>
  auto format(const std::complex<T> &x, FormatCtx &ctx) -> decltype(ctx.out()) {
    format_to(ctx.out(), "(");
    if (style_ == style::pair) {
      base::format(x.real(), ctx);
      format_to(ctx.out(), ",");
      base::format(x.imag(), ctx);
      return format_to(ctx.out(), ")");
    }
    if (x.real() || !x.imag()) base::format(x.real(), ctx);
    if (x.imag()) {
      if (x.real() && x.imag() >= 0 && specs_.sign != sign::plus) format_to(ctx.out(), "+");
      base::format(x.imag(), ctx);
      if (style_ == style::star)
        format_to(ctx.out(), "*i");
      else
        format_to(ctx.out(), "i");
      if (std::is_same<typename std::decay<T>::type, float>::value) format_to(ctx.out(), "f");
      if (std::is_same<typename std::decay<T>::type, long double>::value) format_to(ctx.out(), "l");
    }
    return format_to(ctx.out(), ")");
  }
};
#endif

#endif // CPPDUALS_DUAL
