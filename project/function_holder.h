#pragma once

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
  FuncHolder() = default;
  FuncHolder(T a, T b) : m_a(a), m_b(b) {};
  T a() const {
    return m_a;
  }
  T b() const {
    return m_b;
  }

  virtual T operator()(const T& x) const = 0;
private:
  T m_a = 0;
  T m_b = 0;
};