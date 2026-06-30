#ifndef BRUSSELATOR_HPP
#define BRUSSELATOR_HPP

#include "grid.hpp"
#include <Eigen/Sparse>
#include <Eigen/Core>

// Build the spatial discretization matrix A for a single species with Neumann boundary conditions
Eigen::SparseMatrix<double> buildMatrixA(const GridInfo& grid);

// Build the combined block-diagonal matrix Atilde = [[D1*A, 0], [0, D2*A]]
Eigen::SparseMatrix<double> buildMatrixAtilde(const GridInfo& grid, const Eigen::SparseMatrix<double>& A);

// Evaluate the right-hand side function f(t, C) = Atilde * C + omega(C)
Eigen::VectorXd computeRHS(const GridInfo& grid, const Eigen::SparseMatrix<double>& Atilde, const Eigen::VectorXd& C);

// Evaluate the residual g(C_new) = C_new - C_old - dt * (Atilde * C_new + omega(C_new))
Eigen::VectorXd computeResidual(const GridInfo& grid, const Eigen::SparseMatrix<double>& Atilde, const Eigen::VectorXd& C_new, const Eigen::VectorXd& C_old);

// Assemble the sparse system Jacobian matrix J = I - dt * (Atilde + d_omega/d_C)
Eigen::SparseMatrix<double> buildJacobian(const GridInfo& grid, const Eigen::SparseMatrix<double>& Atilde, const Eigen::VectorXd& C);

#endif // BRUSSELATOR_HPP
