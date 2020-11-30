#pragma once

#include <vector>
#include <sys/types.h>
#include <boost/math/constants/constants.hpp>
#include <cmath>

#include "function_holder.h"
#include "system_solver.h"

namespace bmc = boost::math::constants;

/*
 * synopsis
 */
template <class T>
int generateUniformGrid(ssize_t nodesNumber, T a, T b, std::vector<T>& grid);

template <class T>
int generateChebyshevGrid(ssize_t nodesNumber, T a, T b, std::vector<T>& grid);

///*
// * following classes extends function holder
// */
//template <class T>
//class LagIntCoef;
//
//template <class T>
//class LagInt;
//
//template <class T>
//class CubicSplintIntPart;
//
//template <class T>
//class CubicSplineInt;

/*
 * implementation
 */
template <class T>
int generateUniformGrid(ssize_t nodesNumber, T a, T b, std::vector<T>& grid) {
  grid.resize(nodesNumber + 1);

  if (nodesNumber < 1) {
    return -1;
  }
  T h = (b - a) / static_cast<T>(nodesNumber);

  for (ssize_t i = 0; i < nodesNumber + 1; ++i) {
    grid[i] = a + h * i;
  }

  return 0;
}

template <class T>
int generateChebyshevGrid(ssize_t nodesNumber, T a, T b, std::vector<T>& grid) {
  grid.resize(nodesNumber + 1);

  if (nodesNumber < 1) {
    return -1;
  }

  ssize_t n = nodesNumber;

  for (size_t i = 0; i < n + 1; ++i) {
    grid[i] = (a + b) / 2. + (b - a) / 2. * std::cos((2. * i + 1.) / (2. * (n + 1.)) * bmc::pi<T>());
  }

  return 0;
}

template <class T>
class LagIntCoef : public FuncHolder<T> {
public:
  LagIntCoef() = default;
  LagIntCoef(T a, T b, const std::vector<T>& xi, ssize_t k) : FuncHolder<T>(a, b), m_xi(xi),
                                                              m_n(xi.size()), m_k(k) {};

  T operator()(const T& x) const override {
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
class LagInt : public FuncHolder<T> {
public:
  LagInt(const FuncHolder<T>& func, const std::vector<T>& xi) : FuncHolder<T>(func.a(), func.b()) {
    ssize_t n = xi.size();

    m_ci.resize(n);
    m_yi.resize(n);
    for (ssize_t i = 0; i < n; ++i) {
      m_ci[i] = LagIntCoef<T>(func.a(), func.b(), xi, i);
      m_yi[i] = func(xi[i]);
    }
  };

  T operator()(const T& x) const override {
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
class CubicSplintIntPart : public FuncHolder<T> {
public:
  CubicSplintIntPart(const FuncHolder<T>& func, T a, T b, T c, T d, T xi)
      : FuncHolder<T>(func.a(), func.b()), m_a(a), m_b(b), m_c(c), m_d(d), m_xi(xi) {};

  T operator()(const T& x) const override {
    return m_a + m_b * (x - m_xi) + m_c * std::abs(std::pow(x - m_xi, 2)) + m_d * std::abs(std::pow(x - m_xi, 3));
  }
private:
  T m_a, m_b, m_c, m_d, m_xi;
};

template <class T>
class CubicSplineInt : public FuncHolder<T> {
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

    for (ssize_t i = 0; i < n; ++i) {
      ai[i] = yi(i);
    }

    ub::matrix<T> A; //(ub::zero_matrix(3, n - 2));
    A.resize(3, n - 2);

    A(1, 0) = 2 * (hi(1) + hi(2));
    A(2, 0) = hi(1);
    for (ssize_t i = 1; i < n - 2; ++i) {
      A(0, i) = hi(i + 1);
      A(1, i) = 2 * (hi(i + 1) + hi(i + 2));

      if (i != n - 3) {
        A(2, i) = hi(i + 2);
      }
    }

    ub::matrix<T> B(n - 2, 1);
    for (ssize_t i = 0; i < n - 2; ++i) {
      B(i, 0) = 3 * (gi(i + 2) - gi(i + 1));
    }

    ub::matrix<T> solution;
    tridiagonalMatrixSolve(A, B, solution);

//#define DEBUG_SPLINE_INTERPOLATION
#ifdef DEBUG_SPLINE_INTERPOLATION
    ub::matrix<T> Atest = ub::zero_matrix(n - 2, n - 2);
    for (size_t i = 0; i < n - 2; ++i) {
      if (i != 0) {
        Atest(i, i - 1) = A(0, i);
      }
      Atest(i, i) = A(1, i);
      if (i != n - 3) {
        Atest(i, i + 1) = A(2, i);
      }
    }
    ub::matrix<T> testSolution;
    gaussSolve(Atest, B, testSolution);

    for (ssize_t i = 0; i < testSolution.size1(); ++i) {
      if (std::abs(testSolution(i, 0) - solution(i, 0)) > std::numeric_limits<T>::epsilon()) {
        std::cerr << "warning. tridiagonal method linear equation solution seems to be wrong "
                  << std::endl;
      }
    }
#endif

    ci[n - 1] = ci[0] = 0;
    for (ssize_t k = 1; k < n - 1; ++k) {
      ci[k] = solution(k - 1, 0);
    }

    for (ssize_t i = 0; i < n; ++i) {
      bi[i] = gi(i + 1) - (ci[i + 1] + 2. * ci[i]) * hi(i + 1) / 3.;
      di.at(i) = (ci[i + 1] - ci[i]) / (3. * hi(i + 1));
    }

    for (ssize_t i = 0; i < n - 1; ++i) {
      m_splineInt.emplace_back(func, ai[i], bi[i], ci[i], di[i], xi[i]);
    }

#ifdef DEBUG_SPLINE_INTERPOLATION
    for (ssize_t i = 0; i < n; ++i) {
      std::cout << "a: " << ai[i] << " b: " << bi[i] << " c: " << ci[i] << " d: " << di[i] << std::endl;
    }
#endif
  };

  T operator()(const T& x) const override {
    /*
     * can be implemented using some kind of interval tree
     * in that case worst case complexity will reduced O(n) -> O(log(n))
     */
    for (size_t i = 1; i < m_xi.size(); ++i) {
      if (x < m_xi[i]) {
        return m_splineInt[i - 1](x);
      }
    }

    return m_splineInt.back()(x);
  }

private:
  std::vector<T> m_xi;
  std::vector<CubicSplintIntPart<T>> m_splineInt;
};
