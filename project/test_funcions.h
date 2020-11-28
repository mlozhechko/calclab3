#pragma once

#include <cmath>
#include "function_holder.h"
/*
 * test functions
 */

template <class T>
class Test1 : public FuncHolder<T> {
public:
  Test1() : FuncHolder<T>(-1, 1) {};

  T operator()(const T& x) const override {
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

  T operator()(const T& x) const override {
    return 1. / (1. + static_cast<T>(k) * x * x);
  }
};
template <class T>
using Test2 = Runge<T, 1>;

template <class T>
class Test3 : public FuncHolder<T> {
public:
  Test3() : FuncHolder<T>(-3, 3) {};

  T operator()(const T& x) const override {
    T div = std::abs(std::atan<T>(1. + 10. * x * x));
    return 1. / div;
  }
};