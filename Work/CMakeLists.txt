cmake_minimum_required(VERSION 3.12)
project(BrusselatorSolver CXX)

# Require C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3 package
find_package(Eigen3 REQUIRED)

# Header search paths
include_directories(include)

# Shared source files for the solver and tests
set(SOLVER_SOURCES
    src/brusselator.cpp
    src/integrator.cpp
    src/writer.cpp
)

# Main executable target
add_executable(brusselator src/solver.cpp ${SOLVER_SOURCES})
target_link_libraries(brusselator PRIVATE Eigen3::Eigen)

# Test executable target
add_executable(test_solver tests/test_solver.cpp ${SOLVER_SOURCES})
target_link_libraries(test_solver PRIVATE Eigen3::Eigen)

# Installation rules (installs to CMAKE_INSTALL_PREFIX/bin)
install(TARGETS brusselator test_solver DESTINATION bin)
