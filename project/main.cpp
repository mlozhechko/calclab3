#include <iostream>
#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <sstream>
#include <memory>

#include "interpolation.h"
#include "test_funcions.h"

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
namespace ub = boost::numeric::ublas;

template <class T>
std::vector<T> evalForGrid(const FuncHolder<T>& func, const std::vector<T>& grid, bool isPrint = false) {
  const ssize_t n = grid.size();
  std::vector<T> yValues(n);
  for (size_t i = 0; i < n; ++i) {
    yValues[i] = func(grid[i]);
    if (isPrint) {
      std::cout << grid[i] << " " << func(grid[i]) << std::endl;
    }
  }

  return yValues;
}

template <class T>
int lab3Main(int argc, char** argv) {
  if (!cmdOptionExists(argv, argv + argc, "-grid")) {
    std::cerr << "grid type is not specified" << std::endl;
    return -1;
  }
  std::string gridType = getCmdOption(argv, argv + argc, "-grid");

  if (!cmdOptionExists(argv, argv + argc, "-fid")) {
    std::cerr << "grid type is not specified" << std::endl;
    return -1;
  }
  std::string functionId = getCmdOption(argv, argv + argc, "-fid");
  /*
   * could be implemented better with XXXmap container
   */
  std::shared_ptr<FuncHolder<T>> func{nullptr};
  if (functionId == "test1") {
    func = std::make_shared<Test1<T>>();
  } else if (functionId == "test2") {
    func = std::make_shared<Test2<T>>();
  } else if (functionId == "test3") {
    func = std::make_shared<Test3<T>>();
  } else if (functionId == "runge") {
    func = std::make_shared<Runge<T, 25>>();
  } else if (functionId == "target") {
    func = std::make_shared<TestVar<T>>();
  }

  if (!cmdOptionExists(argv, argv + argc, "-n")) {
    std::cerr << "grid type is not specified" << std::endl;
    return -1;
  }

  ssize_t gridN = std::stoi(getCmdOption(argv, argv + argc, "-n"));

  std::vector<T> grid;
  if (gridType == "uniform") {
    generateUniformGrid(gridN, func->a(), func->b(), grid);
  } else if (gridType == "chebyshev") {
    generateChebyshevGrid(gridN, func->a(), func->b(), grid);
  }

#define DEBUG_RESULTS
#ifdef DEBUG_RESULTS
  std::ostringstream os;
  for (auto& it : grid) {
    os << it << " ";
  }
  std::cout << "generated grid: " << os.str() << std::endl;
#endif

  std::shared_ptr<FuncHolder<T>> interpFunc{nullptr};

  std::string interpolationMethod;
  if (!cmdOptionExists(argv, argv + argc, "-method")) {
    std::cerr << "interpolation method has not been specified" << std::endl;
    return -1;
  }
  interpolationMethod = getCmdOption(argv, argv + argc, "-method");

  if (interpolationMethod == "lagrange") {
    interpFunc = std::make_shared<LagInt<T>>(*func, grid);
  } else if (interpolationMethod == "spline") {
    interpFunc = std::make_shared<CubicSplineInt<T>>(*func, grid);
  } else {
    std::cerr << "interpolation method incorrect" << std::endl;
    return -1;
  }

#ifdef DEBUG_RESULTS
  std::cout << "print souce (x_i, y_i) pairs: " << std::endl;
  evalForGrid(*func, grid, true);
#endif

  std::vector<T> canonicGrid;
  generateUniformGrid(512, func->a(), func->b(), canonicGrid);
  plt::plot(canonicGrid, evalForGrid(*func, canonicGrid), "");
  plt::plot(canonicGrid, evalForGrid(*interpFunc, canonicGrid), "--");
  plt::plot(grid, evalForGrid(*func, grid), "o");

  plt::ylim(-1, 2);

//  std::string filename = "./results/" + functionId + "_" + gridType + "_" + std::to_string(gridN) + ".png";
  std::string filename = "results/result.png";

  std::cout << "calculation error norm: "
            << calculationErrorNorm(*func, *interpFunc, canonicGrid)
            << std::endl;

  plt::save(filename, 500);

  return 0;
}


int main(int argc, char** argv) {
  if (!cmdOptionExists(argv, argv + argc, "-precision")) {
    std::cerr << "precision is not specified" << std::endl;
    return -1;
  }

  std::string precision = getCmdOption(argv, argv + argc, "-precision");
  if (precision == "double") {
    return lab3Main<double>(argc, argv);
  } else if (precision == "float") {
    return lab3Main<float>(argc, argv);
  }

  std::cerr << "precision cannot be parsed correctly" << std::endl;
  return -2;
}