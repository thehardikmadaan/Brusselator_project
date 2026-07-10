#include "integrator.hpp"
#include "brusselator.hpp"
#include <iostream>
#include <Eigen/Dense>

void stepRK4(const GridInfo& grid, const Eigen::SparseMatrix<double>& Atilde, Eigen::VectorXd& C) {
    double dt = grid.dt;

    Eigen::VectorXd k1 = computeRHS(grid, Atilde, C);
    Eigen::VectorXd k2 = computeRHS(grid, Atilde, C + 0.5 * dt * k1);
    Eigen::VectorXd k3 = computeRHS(grid, Atilde, C + 0.5 * dt * k2);
    Eigen::VectorXd k4 = computeRHS(grid, Atilde, C + dt * k3);

    C += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

void stepExplicitEuler(const GridInfo& grid, const Eigen::SparseMatrix<double>& Atilde, Eigen::VectorXd& C) {
    C += grid.dt * computeRHS(grid, Atilde, C);
}

void stepImplicitEuler(const GridInfo& grid, const Eigen::SparseMatrix<double>& Atilde, Eigen::VectorXd& C, int& total_newton_its) {
    double tol = 1e-6;
    int max_its = 50;
    
    Eigen::VectorXd C_old = C;
    Eigen::VectorXd C_new = C;
    
    // Build and factorize the Jacobian once at the start of the time step
    Eigen::SparseMatrix<double> J = buildJacobian(grid, Atilde, C_new);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(J);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Jacobian factorization failed!\n";
        return;
    }
    
    int it = 0;
    bool converged = false;
    
    while (it < max_its) {
        Eigen::VectorXd g = computeResidual(grid, Atilde, C_new, C_old);
        double residual_norm = g.lpNorm<Eigen::Infinity>();
        
        if (residual_norm < tol) {
            converged = true;
            break;
        }
        
        Eigen::VectorXd delta_C = solver.solve(-g);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Linear solve failed in Newton's method!\n";
            return;
        }
        
        C_new += delta_C;
        it++;
    }
    
    total_newton_its += it;
    
    if (!converged) {
        std::cerr << "Newton's method failed to converge! Residual norm: " 
                  << computeResidual(grid, Atilde, C_new, C_old).lpNorm<Eigen::Infinity>() << "\n";
    }
    
    C = C_new;
}
