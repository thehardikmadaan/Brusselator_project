cmake_minimum_required(VERSION 3.12)
project(BrusselatorSolver CXX)

# Require C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3 package
find_package(Eigen3 REQUIRED)

# Find VTK package
find_package(VTK REQUIRED)

# Header search paths
include_directories(include)

# Shared source files for the solver and tests
set(SOLVER_SOURCES
    src/brusselator.cpp
    src/integrator.cpp
    src/writer.cpp
)

# Create a library for the core logic
add_library(brusselator_lib STATIC ${SOLVER_SOURCES})
target_link_libraries(brusselator_lib PUBLIC Eigen3::Eigen ${VTK_LIBRARIES})

# Main executable target
add_executable(brusselator src/solver.cpp)
target_link_libraries(brusselator PRIVATE brusselator_lib)

# Test executable target
add_executable(test_solver tests/test_solver.cpp)
target_link_libraries(test_solver PRIVATE brusselator_lib)

# Installation rules (installs to CMAKE_INSTALL_PREFIX/bin)
install(TARGETS brusselator test_solver DESTINATION bin)
