cmake_minimum_required(VERSION 3.10)
project(calc4)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Python3 COMPONENTS Development NumPy)

include_directories(project)
add_executable(calc4 project/main.cpp project/system_solver.h project/utils.h
        project/matrix_utils.h project/eigen_solver.h project/stat_holder.h
        project/interpolation.h project/function_holder.h project/test_funcions.h)

target_include_directories(calc4 PRIVATE ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(calc4 Python3::Python Python3::NumPy)