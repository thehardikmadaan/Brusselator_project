#include "grid.hpp"
#include "brusselator.hpp"
#include "integrator.hpp"
#include <iostream>
#include <cmath>
#include <cassert>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vtkSmartPointer.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>

// Test 1: Verify that the row sums of the Laplacian matrix A are exactly 0 (Neumann BCs)
void testMatrixCorrectness() {
    std::cout << "Running Test 1: Matrix Correctness (Row Sums)... ";
    
    GridInfo grid;
    grid.nx = 10;
    grid.ny = 10;
    
    Eigen::SparseMatrix<double> A = buildMatrixA(grid);
    int N = grid.nx * grid.ny;
    
    for (int i = 0; i < N; ++i) {
        double row_sum = 0.0;
        for (int k = 0; k < A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
                if (it.row() == i) {
                    row_sum += it.value();
                }
            }
        }
        // Row sum should be 0 up to machine precision
        if (std::abs(row_sum) > 1e-11) {
            std::cerr << "\nFAIL: Row " << i << " sum is " << row_sum << " (expected 0.0)\n";
            std::exit(1);
        }
    }
    
    std::cout << "PASSED\n";
}

// Test 2: Verify mass conservation when reactions are disabled
void testMassConservation() {
    std::cout << "Running Test 2: Mass Conservation (Reactions Off)... ";
    
    GridInfo grid;
    grid.nx = 20;
    grid.ny = 20;
    grid.duration = 0.1;
    grid.dt = 0.01;
    grid.reactions = false; // Disable reactions
    
    int N = grid.nx * grid.ny;
    Eigen::VectorXd C1 = Eigen::VectorXd::Constant(N, 1.0);
    Eigen::VectorXd C2 = Eigen::VectorXd::Constant(N, 2.0);
    
    // Add some random perturbations to make it interesting
    C1[5 * grid.nx + 5] = 1.5;
    C2[15 * grid.nx + 15] = 2.5;
    
    double initial_mass1 = C1.sum();
    double initial_mass2 = C2.sum();
    
    Eigen::VectorXd C(2 * N);
    C.head(N) = C1;
    C.tail(N) = C2;
    
    Eigen::SparseMatrix<double> A = buildMatrixA(grid);
    Eigen::SparseMatrix<double> Atilde = buildMatrixAtilde(grid, A);
    
    // Run 10 steps of RK4
    double t = 0.0;
    while (t < grid.duration) {
        stepRK4(grid, Atilde, C);
        t += grid.dt;
    }
    
    double final_mass1 = C.head(N).sum();
    double final_mass2 = C.tail(N).sum();
    
    // Mass should be conserved exactly (up to tiny rounding errors)
    if (std::abs(final_mass1 - initial_mass1) > 1e-10 || std::abs(final_mass2 - initial_mass2) > 1e-10) {
        std::cerr << "\nFAIL: Mass not conserved!\n"
                  << "  C1 Initial: " << initial_mass1 << " | Final: " << final_mass1 << "\n"
                  << "  C2 Initial: " << initial_mass2 << " | Final: " << final_mass2 << "\n";
        std::exit(1);
    }
    
    std::cout << "PASSED\n";
}

// Test 3: Verify that RK4 and Implicit Euler agree for a small time step
void testSolverAgreement() {
    std::cout << "Running Test 3: Solver Agreement (RK4 vs Implicit Euler)... ";
    
    GridInfo grid;
    grid.nx = 15;
    grid.ny = 15;
    grid.dt = 0.001; // very small time step
    grid.Ca = 1.0;
    grid.Cb = 3.0;
    grid.reactions = true;
    
    int N = grid.nx * grid.ny;
    Eigen::VectorXd C_init(2 * N);
    C_init.head(N) = Eigen::VectorXd::Constant(N, grid.Ca);
    C_init.tail(N) = Eigen::VectorXd::Constant(N, grid.Cb / grid.Ca);
    
    // Add perturbations
    C_init[3 * grid.nx + 3] *= 1.2;
    C_init[12 * grid.nx + 12] *= 0.9;
    
    Eigen::SparseMatrix<double> A = buildMatrixA(grid);
    Eigen::SparseMatrix<double> Atilde = buildMatrixAtilde(grid, A);
    
    // 1. Take one step of RK4
    Eigen::VectorXd C_rk4 = C_init;
    stepRK4(grid, Atilde, C_rk4);
    
    // 2. Take one step of Implicit Euler
    Eigen::VectorXd C_implicit = C_init;
    int newton_its = 0;
    stepImplicitEuler(grid, Atilde, C_implicit, newton_its);
    
    // Check difference between solutions
    double diff = (C_rk4 - C_implicit).lpNorm<Eigen::Infinity>();
    
    // For dt = 0.001, RK4 (O(dt^4)) and Implicit Euler (O(dt)) should be extremely close
    if (diff > 1e-5) {
        std::cerr << "\nFAIL: Solvers disagree! Max difference: " << diff << " (expected < 1e-5)\n";
        std::exit(1);
    }
    
    std::cout << "PASSED\n";
}

int main() {
    std::cout << "==================================================\n";
    std::cout << "Running Brusselator Solver Unit Tests...\n";
    std::cout << "==================================================\n";
    
    testMatrixCorrectness();
    testMassConservation();
    testSolverAgreement();
    
    std::cout << "==================================================\n";
    std::cout << "All unit tests completed successfully!\n";
    std::cout << "==================================================\n";
    return 0;
}
