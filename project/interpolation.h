#pragma once

#include <vector>
#include <sys/types.h>
#include <boost/math/constants/constants.hpp>
#include <cmath>

#include "function_holder.h"

namespace bmc = boost::math::constants;

/*
 * synopsis
 */
template <class T>
int generateUniformGrid(ssize_t nodesNumber, T a, T b, std::vector<T>& grid);

template <class T>
int generateChebyshevGrid(ssize_t nodesNumber, T a, T b, std::vector<T>& grid);

/*
 * implementation
 */
template <class T>
int generateUniformGrid(ssize_t nodesNumber, T a, T b, std::vector<T>& grid) {
  grid.resize(nodesNumber);

  if (nodesNumber < 1) {
    return -1;
  }
  T h = (b - a) / static_cast<T>(nodesNumber);

  for (ssize_t i = 0; i < nodesNumber; ++i) {
    grid[i] = a + h * i;
  }

  return 0;
}

template <class T>
int generateChebyshevGrid(ssize_t nodesNumber, T a, T b, std::vector<T>& grid) {
  grid.resize(nodesNumber);

  if (nodesNumber < 1) {
    return -1;
  }

  ssize_t n = nodesNumber;

  for (size_t i = 0; i < n; ++i) {
    grid[i] = (a + b) / 2. + (b - a) / 2. + std::cos((2. * i + 1.) / (2. * (n + 1.)) * bmc::pi<T>());
  }

  return 0;
}

