#pragma once

#include <cmath>
#include "system_solver.h"
#include "boost/icl/interval_map.hpp"

namespace icl = boost::icl;
/*
 * we'll use polymorphism to easily pass functions as
 * arguments. it's not the best but the fasted way to implement that
 * kind of work
 */

/*
 * function base
 */

template <class T>
class FuncHolder {
public:
  FuncHolder(T a, T b) : m_a(a), m_b(b) {};
  T a() const {
    return m_a;
  }
  T b() const {
    return m_b;
  }

  virtual T operator()(const T& x) = 0;
private:
  T m_a = 0;
  T m_b = 0;
};

/*
 * test functions
 */

template <class T>
class Test1 : public FuncHolder<T> {
public:
  Test1() : FuncHolder<T>(-1, 1) {};

  T operator()(const T& x) override {
    return x * x;
  }
};

/*
 * Runge's function
 */
template <class T, int k>
class Runge : public FuncHolder<T> {
public:
  Runge() : FuncHolder<T>(-1, 1) {};

  T operator()(const T& x) override {
    return 1. / (1. + static_cast<T>(k) * x * x);
  }
};
template <class T>
using Test2 = Runge<T, 1>;

template <class T>
class Test3 : public FuncHolder<T> {
public:
  Test3() : FuncHolder<T>(-3, 3) {};

  T operator()(const T& x) override {
    T div = std::abs(std::atan<T>(1. + 10. * x * x));
    return 1. / div;
  }
};

template <class T>
class LagIntCoef : public FuncHolder<T> {
public:
  LagIntCoef(T a, T b, const std::vector<T>& xi, ssize_t k) : FuncHolder<T>(a, b), m_xi(xi),
                                                              m_n(xi.size()), m_k(k) {};

  T operator()(const T& x) override {
    T acc = 1;

    for (ssize_t i = 0; i < m_n; ++i) {
      if (i == m_k) {
        continue;
      }

      acc *= (x - m_xi[i]) / (m_xi[m_k] - m_xi[i]);
    }

    return acc;
  }
private:
  std::vector<T> m_xi{};
  ssize_t m_n = 0;
  ssize_t m_k = 0;
};

template <class T>
class LagInt : FuncHolder<T> {
public:
  LagInt(const FuncHolder<T>& func, const std::vector<T>& xi) : FuncHolder<T>(func.a(), func.b()) {
    ssize_t n = xi.size();

    m_ci.resize(n);
    m_yi.resize(n);
    for (ssize_t i = 0; i < n; ++i) {
      m_ci = LagIntCoef<T>(func.a(), func.b(), xi, i);
      m_yi = func(xi[i]);
    }
  };

  T operator()(const T& x) {
    ssize_t n = m_yi.size();
    T acc = 0;
    for (ssize_t i = 0; i < n; ++i) {
      acc += m_ci[i](x) * m_yi[i];
    }

    return acc;
  }

private:
  std::vector<LagIntCoef<T>> m_ci;
  std::vector<T> m_yi;
};

template <class T>
class CubicSplintIntPart : FuncHolder<T> {
public:
  CubicSplintIntPart(const FuncHolder<T>& func, T a, T b, T c, T d, T xi)
  : FuncHolder<T>(func.a(), func.b()), m_a(a), m_b(b), m_c(c), m_d(d), m_xi(xi) {};

  T operator()(const T& x) {
    return m_a + m_b * (x - m_xi) + m_c * std::abs(std::pow(x - m_xi, 2)) + m_d * std::abs(std::pow(x - m_xi, 3));
  }
private:
  T m_a, m_b, m_c, m_d, m_xi;
};

template <class T>
class CubicSplineInt : FuncHolder<T> {
public:
  CubicSplineInt(const FuncHolder<T>& func, const std::vector<T> xi)
      : FuncHolder<T>(func.a(), func.b()) {
    m_xi = xi;
    const ssize_t n = xi.size();

    std::vector<T> ai(n);
    std::vector<T> bi(n);
    std::vector<T> ci(n);
    std::vector<T> di(n);

    auto hi = [&](ssize_t i) -> T {
      return xi[i] - xi[i - 1];
    };
    auto yi = [&](ssize_t i) -> T {
      return func(xi[i]);
    };
    auto gi = [&](ssize_t i) -> T {
      return (yi(i) - yi(i - 1)) / hi(i);
    };

    for (ssize_t i = 1; i < n; ++i) {
      ai[i] = yi(i);
    }

    ub::matrix<T> A = ub::zero_matrix(3, n - 2);

    A(1, 0) = 2 * (hi[1] + hi[2]);
    A(2, 0) = hi[1];
    for (ssize_t i = 1; i < n - 2; ++i) {
      A(0, i) = hi(i + 1);
      A(1, i) = 2 * (hi(i + 1) + hi(i + 2));
      A(2, i) = hi(i + 2);
    }

    ub::matrix<T> B(n - 2, 1);
    for (ssize_t i = 0; i < n - 2; ++i) {
      B(i, 0) = 3 * (gi(i + 2) - gi(i + 1));
    }

    ub::matrix<T> solution;
    tridiagonalMatrixSolve(A, B, solution);

    ci[n] = ci[0] = 0;
    for (ssize_t k = 1; k < n - 1; ++k) {
      ci[k] = solution(k - 1, 0);
    }

    for (ssize_t i = 0; i < n; ++i) {
      bi[i] = gi(i + 1) - (ci[i + 1] + 2. * ci[i]) * hi(i + 1) / 3.;
      di[i] = (ci[i + 1] - ci[i]) / (3. * hi(i + 1));
    }

    for (ssize_t i = 0; i < n - 1; ++i) {
      m_splineInt.emplace_back(func, ai[i], bi[i], di[i], ci[i], xi[i]);
    }
  };

  T operator()(const T& x) {
    /*
     * can be implemented using some kind of interval tree
     * in that case worst case complexity will reduced O(n) -> O(log(n))
     */
    for (size_t i = 1; i < m_xi.size(); ++i) {
      if (x < m_xi[i]) {
        return m_splineInt[i - 1](x);
      }
    }
  }

private:
  std::vector<T> m_xi;
  std::vector<CubicSplintIntPart<T>> m_splineInt;
};