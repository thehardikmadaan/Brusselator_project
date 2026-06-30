cmake_minimum_required(VERSION 3.12)
project(BrusselatorSolver CXX)

# Require C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3 package
find_package(Eigen3 REQUIRED)

# Add the executable
add_executable(brusselator main.cpp)

# Link the required libraries (Eigen3)
target_link_libraries(brusselator PRIVATE Eigen3::Eigen)

# Installation rules
install(TARGETS brusselator DESTINATION bin)
