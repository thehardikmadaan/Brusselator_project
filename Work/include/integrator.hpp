#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include "grid.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>

// Perform a single explicit RK4 time step to advance C in place
void stepRK4(const GridInfo& grid, const Eigen::SparseMatrix<double>& Atilde, Eigen::VectorXd& C);

// Perform a single implicit Euler time step using Newton's method
void stepImplicitEuler(const GridInfo& grid, const Eigen::SparseMatrix<double>& Atilde, Eigen::VectorXd& C, int& total_newton_its);

#endif // INTEGRATOR_HPP
