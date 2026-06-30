#include <iostream>
#include <cstdlib>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "brusselator.h"

// Assemble the spatial discretization matrix A for a single species
Eigen::SparseMatrix<double> buildMatrixA(const SimParams& params) {
    int N = params.nx * params.ny;
    Eigen::SparseMatrix<double> A(N, N);
    std::vector<Eigen::Triplet<double>> triplets;

    double dx_sq = params.dx * params.dx;
    double factor = 1.0 / dx_sq;

    // Loop over the entire grid
    for (int i = 0; i < params.ny; ++i) {
        for (int j = 0; j < params.nx; ++j) {
            int k = i * params.nx + j;

            // Center point: -4 / (dx^2)
            triplets.push_back(Eigen::Triplet<double>(k, k, -4.0 * factor));

            // West (j - 1)
            if (j > 0) {
                triplets.push_back(Eigen::Triplet<double>(k, k - 1, factor));
            } else {
                triplets.push_back(Eigen::Triplet<double>(k, k + 1, factor)); // Neumann BC reflection
            }

            // East (j + 1)
            if (j < params.nx - 1) {
                triplets.push_back(Eigen::Triplet<double>(k, k + 1, factor));
            } else {
                triplets.push_back(Eigen::Triplet<double>(k, k - 1, factor)); // Neumann BC reflection
            }

            // South (i - 1)
            if (i > 0) {
                triplets.push_back(Eigen::Triplet<double>(k, k - params.nx, factor));
            } else {
                triplets.push_back(Eigen::Triplet<double>(k, k + params.nx, factor)); // Neumann BC reflection
            }

            // North (i + 1)
            if (i < params.ny - 1) {
                triplets.push_back(Eigen::Triplet<double>(k, k + params.nx, factor));
            } else {
                triplets.push_back(Eigen::Triplet<double>(k, k - params.nx, factor)); // Neumann BC reflection
            }
        }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

// Assemble the combined block-diagonal matrix Atilde
Eigen::SparseMatrix<double> buildMatrixAtilde(const SimParams& params, const Eigen::SparseMatrix<double>& A) {
    int N = params.nx * params.ny;
    Eigen::SparseMatrix<double> Atilde(2 * N, 2 * N);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(A.nonZeros() * 2);

    // Loop over the non-zero elements of A to build the diagonal blocks
    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            // Top-left block: D1 * A (activator)
            triplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), params.D1 * it.value()));
            // Bottom-right block: D2 * A (inhibitor)
            triplets.push_back(Eigen::Triplet<double>(it.row() + N, it.col() + N, params.D2 * it.value()));
        }
    }

    Atilde.setFromTriplets(triplets.begin(), triplets.end());
    return Atilde;
}

// Evaluate the RHS function f(t, C) = Atilde * C + omega(C)
Eigen::VectorXd computeRHS(const SimParams& params, const Eigen::SparseMatrix<double>& Atilde, const Eigen::VectorXd& C) {
    int N = params.nx * params.ny;
    
    // 1. Compute diffusion term: Atilde * C
    Eigen::VectorXd rhs = Atilde * C;

    // 2. Add nonlinear reaction terms element-wise
    for (int i = 0; i < N; ++i) {
        double c1 = C[i];          // Activator concentration at point i
        double c2 = C[i + N];      // Inhibitor concentration at point i
        double c1_sq = c1 * c1;

        // omega_1 = Ca + c1^2 * c2 - Cb * c1 - c1
        double omega1 = params.Ca + c1_sq * c2 - params.Cb * c1 - c1;
        // omega_2 = Cb * c1 - c1^2 * c2
        double omega2 = params.Cb * c1 - c1_sq * c2;

        rhs[i] += omega1;
        rhs[i + N] += omega2;
    }

    return rhs;
}

// Perform a single explicit RK4 time step
void stepRK4(const SimParams& params, const Eigen::SparseMatrix<double>& Atilde, Eigen::VectorXd& C) {
    double dt = params.dt;

    Eigen::VectorXd k1 = computeRHS(params, Atilde, C);
    Eigen::VectorXd k2 = computeRHS(params, Atilde, C + 0.5 * dt * k1);
    Eigen::VectorXd k3 = computeRHS(params, Atilde, C + 0.5 * dt * k2);
    Eigen::VectorXd k4 = computeRHS(params, Atilde, C + dt * k3);

    C += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

// Evaluate the residual g(C_new) = C_new - C_old - dt * (Atilde * C_new + omega(C_new))
Eigen::VectorXd computeResidual(const SimParams& params, const Eigen::SparseMatrix<double>& Atilde, const Eigen::VectorXd& C_new, const Eigen::VectorXd& C_old) {
    int N = params.nx * params.ny;
    
    // g = C_new - C_old - dt * (Atilde * C_new)
    Eigen::VectorXd g = C_new - C_old - params.dt * (Atilde * C_new);
    
    // Subtract reaction term contribution: - dt * omega(C_new)
    for (int i = 0; i < N; ++i) {
        double c1 = C_new[i];
        double c2 = C_new[i + N];
        double c1_sq = c1 * c1;
        
        double omega1 = params.Ca + c1_sq * c2 - params.Cb * c1 - c1;
        double omega2 = params.Cb * c1 - c1_sq * c2;
        
        g[i] -= params.dt * omega1;
        g[i + N] -= params.dt * omega2;
    }
    
    return g;
}

// Assemble the sparse Jacobian matrix J = I - dt * (Atilde + d_omega/d_C)
Eigen::SparseMatrix<double> buildJacobian(const SimParams& params, const Eigen::SparseMatrix<double>& Atilde, const Eigen::VectorXd& C) {
    int N = params.nx * params.ny;
    Eigen::SparseMatrix<double> J(2 * N, 2 * N);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(Atilde.nonZeros() + 4 * N);
    
    // Add contribution of: I - dt * Atilde
    for (int k = 0; k < Atilde.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Atilde, k); it; ++it) {
            int r = it.row();
            int c = it.col();
            double val = -params.dt * it.value();
            if (r == c) {
                val += 1.0; // Add Identity diagonal
            }
            triplets.push_back(Eigen::Triplet<double>(r, c, val));
        }
    }
    
    // Add local reaction derivative contributions: -dt * d_omega/d_C
    for (int p = 0; p < N; ++p) {
        double c1 = C[p];
        double c2 = C[p + N];
        
        double d_w1_dc1 = 2.0 * c1 * c2 - params.Cb - 1.0;
        double d_w1_dc2 = c1 * c1;
        double d_w2_dc1 = params.Cb - 2.0 * c1 * c2;
        double d_w2_dc2 = -c1 * c1;
        
        triplets.push_back(Eigen::Triplet<double>(p, p, -params.dt * d_w1_dc1));
        triplets.push_back(Eigen::Triplet<double>(p, p + N, -params.dt * d_w1_dc2));
        triplets.push_back(Eigen::Triplet<double>(p + N, p, -params.dt * d_w2_dc1));
        triplets.push_back(Eigen::Triplet<double>(p + N, p + N, -params.dt * d_w2_dc2));
    }
    
    J.setFromTriplets(triplets.begin(), triplets.end());
    return J;
}

// Perform a single implicit Euler step using Inexact Newton's method
void stepImplicitEuler(const SimParams& params, const Eigen::SparseMatrix<double>& Atilde, Eigen::VectorXd& C, int& total_newton_its) {
    double tol = 1e-6;
    int max_its = 50;
    
    Eigen::VectorXd C_old = C; // Solution at time level k
    Eigen::VectorXd C_new = C; // Initial guess for time level k+1 (using C^k)
    
    // 1. Build and decompose the Jacobian ONCE at the start of the time step
    Eigen::SparseMatrix<double> J = buildJacobian(params, Atilde, C_new);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(J);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Jacobian factorization failed!\n";
        return;
    }
    
    int it = 0;
    bool converged = false;
    
    while (it < max_its) {
        // 2. Compute residual
        Eigen::VectorXd g = computeResidual(params, Atilde, C_new, C_old);
        
        // Check convergence using the Infinity norm
        double residual_norm = g.lpNorm<Eigen::Infinity>();
        if (residual_norm < tol) {
            converged = true;
            break;
        }
        
        // 3. Solve the linear system J * delta_C = -g
        Eigen::VectorXd delta_C = solver.solve(-g);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Linear solve failed in Newton's method!\n";
            return;
        }
        
        // 4. Update the solution
        C_new += delta_C;
        it++;
    }
    
    total_newton_its += it;
    
    if (!converged) {
        std::cerr << "Newton's method failed to converge! Residual norm: " 
                  << computeResidual(params, Atilde, C_new, C_old).lpNorm<Eigen::Infinity>() << "\n";
    }
    
    C = C_new;
}

int main(int argc, char *argv[]) {
    // We allow an optional 8th argument to specify the solver method ("rk4" or "implicit")
    std::string method = "rk4";
    if (argc == 9) {
        method = argv[8];
    } else if (argc != 8) {
        std::cerr << "Usage: " << argv[0]
                  << " <nx> <ny> <duration> <dt> <write_interval> <Ca> <Cb> [<method: rk4|implicit>]\n";
        return 1;
    }

    // Parse command line arguments
    SimParams params;
    params.nx = std::atoi(argv[1]);
    params.ny = std::atoi(argv[2]);
    params.duration = std::atof(argv[3]);
    params.dt = std::atof(argv[4]);
    params.write_interval = std::atof(argv[5]);
    params.Ca = std::atof(argv[6]);
    params.Cb = std::atof(argv[7]);

    // Hardcoded physical constants as per project sheet
    params.dx = 0.005;
    params.dy = 0.005;
    params.D1 = 1e-5;
    params.D2 = 1e-6;

    std::cout << "Starting Brusselator Simulation...\n";
    std::cout << "Method: " << (method == "implicit" ? "Implicit Euler (Inexact Newton)" : "Explicit RK4") << "\n";
    std::cout << "Grid: " << params.nx << "x" << params.ny << "\n";
    std::cout << "Domain: " << (params.nx - 1) * params.dx << "m x " << (params.ny - 1) * params.dy << "m\n";

    // Initialize data structures
    int num_points = params.nx * params.ny;
    Eigen::VectorXd C1 = Eigen::VectorXd::Constant(num_points, params.Ca);
    Eigen::VectorXd C2 = Eigen::VectorXd::Constant(num_points, params.Cb / params.Ca);

    // Apply initial perturbations (if grid size allows)
    if (params.nx > 30 && params.ny > 30) {
        int p1_idx = 10 * params.nx + 10;
        int p2_idx = 30 * params.nx + 30;

        C1[p1_idx] *= 1.2;
        C2[p1_idx] *= 1.2;

        C1[p2_idx] *= 0.9;
        C2[p2_idx] *= 0.9;
        std::cout << "Applied initial perturbations at (10,10) and (30,30).\n";
    }

    // Build spatial discretization matrices
    std::cout << "Building sparse matrix A...\n";
    Eigen::SparseMatrix<double> A = buildMatrixA(params);

    std::cout << "Building combined matrix Atilde...\n";
    Eigen::SparseMatrix<double> Atilde = buildMatrixAtilde(params, A);

    // Combine concentrations into a single state vector C
    Eigen::VectorXd C(2 * num_points);
    C.head(num_points) = C1;
    C.tail(num_points) = C2;

    // Time integration loop
    double t = 0.0;
    int step = 0;
    int total_newton_its = 0;

    std::cout << "Starting time integration loop...\n";
    while (t < params.duration) {
        // Adjust time step at the end to match duration exactly
        double current_dt = params.dt;
        if (t + current_dt > params.duration) {
            current_dt = params.duration - t;
        }

        // Advance one time step
        if (method == "implicit") {
            int previous_its = total_newton_its;
            stepImplicitEuler(params, Atilde, C, total_newton_its);
            int step_its = total_newton_its - previous_its;

            t += current_dt;
            step++;

            // Write current time and Newton iterations to stdout
            if (step % 10 == 0 || t >= params.duration) {
                std::cout << "Time: " << t << " s | Step: " << step 
                          << " | Newton Iterations: " << step_its << "\n";
            }
        } else {
            stepRK4(params, Atilde, C);
            t += current_dt;
            step++;

            if (step % 100 == 0 || t >= params.duration) {
                std::cout << "Time: " << t << " s | Step: " << step << "\n";
            }
        }
    }

    std::cout << "Brusselator Simulation completed successfully.\n";
    if (method == "implicit") {
        std::cout << "Total Newton Iterations: " << total_newton_its << "\n";
    }
    return 0;
}
