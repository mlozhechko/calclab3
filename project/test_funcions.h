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

template <class T>
class TestVar : public FuncHolder<T> {
public:
  TestVar() : FuncHolder<T>(-1, 1) {};

  T operator()(const T& x) const override {
//    R1 := Sin((2*Degree(x,2)-x+2*Degree(7,1/3)-5)/2);
//    R2 := Exp((Degree(x,2)+2*x+1)/(7*x+1));

    T r1 = std::sin((2. * std::pow(x, 2.) - x + 2. * std::pow(7., 1. / 3.) - 5.) / 2.);
    T r2 = std::exp((std::pow(x, 2.) + 2. * x + 1.) / (7. * x + 1.));

    return r1 + r2 - 1.5;
  }
};

template <class T>
class ConstTest : public FuncHolder<T> {
public:
  ConstTest() : FuncHolder<T>(-1, 1) {};

  T operator()(const T& x) const override {
    return 1;
  }
};

template <class T>
class LinearTest : public FuncHolder<T> {
public:
  LinearTest() : FuncHolder<T>(-1, 1) {};

  T operator()(const T& x) const override {
    return x;
  }
};

template <class T>
class QuadTest : public FuncHolder<T> {
public:
  QuadTest() : FuncHolder<T>(-1, 1) {};

  T operator()(const T& x) const override {
    return x * x;
  }
};

template <class T>
class SinTest : public FuncHolder<T> {
public:
  SinTest(T a, T b) : FuncHolder<T>(a, b) {};

  T operator()(const T& x) const override {
    return std::sin(bmc::pi<T>() * x);
  }
};