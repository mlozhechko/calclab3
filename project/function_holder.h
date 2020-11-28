#pragma once

#include <cmath>

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
//
//template <class T>
//class CubicSplintIntPart : FuncHolder<T> {
//public:
//  CubicSplintIntPart(const FuncHolder<T>& func, const std::vector<T> xi, ssize_t k)
//  : FuncHolder<T>(func.a(), func.b()) {
//
//  }
//private:
//  T a;
//  T b;
//  T c;
//  T d;
//};

template <class T>
class CubicSplineInt : FuncHolder<T> {
public:
  CubicSplineInt(const FuncHolder<T>& func, const std::vector<T> xi)
      : FuncHolder<T>(func.a(), func.b()) {
    const ssize_t n = xi.size();

    std::vector<T> ai(n);

    for (ssize_t i = 1; i < n; ++i) {
      ai[i] = func(xi[i - 1]);
    }

    // to find c coefficients fast enough we have to use tridiagonal matrix algorithm (see Thomas algorithm)


  };

private:

};