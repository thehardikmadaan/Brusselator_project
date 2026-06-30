#ifndef BRUSSELATOR_H
#define BRUSSELATOR_H

#include <Eigen/Sparse>

// Struct to store all simulation and physical parameters
struct SimParams {
    int nx;                 // Number of grid points in x-direction
    int ny;                 // Number of grid points in y-direction
    double dx;              // Grid spacing in x-direction [m]
    double dy;              // Grid spacing in y-direction [m]
    double duration;        // Total simulation time [s]
    double dt;              // Time step size [s]
    double write_interval;  // Time interval between writing output files [s]
    double Ca;              // Brusselator parameter Ca (reaction rate/concentration)
    double Cb;              // Brusselator parameter Cb (reaction rate/concentration)
    double D1;              // Diffusion coefficient for species 1 [m^2/s]
    double D2;              // Diffusion coefficient for species 2 [m^2/s]
};

// Function declarations
Eigen::SparseMatrix<double> buildMatrixA(const SimParams& params);
Eigen::SparseMatrix<double> buildMatrixAtilde(const SimParams& params, const Eigen::SparseMatrix<double>& A);
Eigen::VectorXd computeRHS(const SimParams& params, const Eigen::SparseMatrix<double>& Atilde, const Eigen::VectorXd& C);
void stepRK4(const SimParams& params, const Eigen::SparseMatrix<double>& Atilde, Eigen::VectorXd& C);
Eigen::VectorXd computeResidual(const SimParams& params, const Eigen::SparseMatrix<double>& Atilde, const Eigen::VectorXd& C_new, const Eigen::VectorXd& C_old);
Eigen::SparseMatrix<double> buildJacobian(const SimParams& params, const Eigen::SparseMatrix<double>& Atilde, const Eigen::VectorXd& C);
void stepImplicitEuler(const SimParams& params, const Eigen::SparseMatrix<double>& Atilde, Eigen::VectorXd& C, int& total_newton_its);

#endif // BRUSSELATOR_H
