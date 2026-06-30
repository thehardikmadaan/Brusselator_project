#include "brusselator.hpp"
#include <vector>

Eigen::SparseMatrix<double> buildMatrixA(const GridInfo& grid) {
    int N = grid.nx * grid.ny;
    Eigen::SparseMatrix<double> A(N, N);
    std::vector<Eigen::Triplet<double>> triplets;

    double dx_sq = grid.dx * grid.dx;
    double factor = 1.0 / dx_sq;

    for (int i = 0; i < grid.ny; ++i) {
        for (int j = 0; j < grid.nx; ++j) {
            int k = i * grid.nx + j;

            int active_neighbors = 0;

            // West (j - 1)
            if (j > 0) {
                triplets.push_back(Eigen::Triplet<double>(k, k - 1, factor));
                active_neighbors++;
            }

            // East (j + 1)
            if (j < grid.nx - 1) {
                triplets.push_back(Eigen::Triplet<double>(k, k + 1, factor));
                active_neighbors++;
            }

            // South (i - 1)
            if (i > 0) {
                triplets.push_back(Eigen::Triplet<double>(k, k - grid.nx, factor));
                active_neighbors++;
            }

            // North (i + 1)
            if (i < grid.ny - 1) {
                triplets.push_back(Eigen::Triplet<double>(k, k + grid.nx, factor));
                active_neighbors++;
            }

            // Center point: -active_neighbors / (dx^2)
            triplets.push_back(Eigen::Triplet<double>(k, k, -static_cast<double>(active_neighbors) * factor));
        }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

Eigen::SparseMatrix<double> buildMatrixAtilde(const GridInfo& grid, const Eigen::SparseMatrix<double>& A) {
    int N = grid.nx * grid.ny;
    Eigen::SparseMatrix<double> Atilde(2 * N, 2 * N);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(A.nonZeros() * 2);

    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            // Top-left block: D1 * A
            triplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), grid.D1 * it.value()));
            // Bottom-right block: D2 * A
            triplets.push_back(Eigen::Triplet<double>(it.row() + N, it.col() + N, grid.D2 * it.value()));
        }
    }

    Atilde.setFromTriplets(triplets.begin(), triplets.end());
    return Atilde;
}

Eigen::VectorXd computeRHS(const GridInfo& grid, const Eigen::SparseMatrix<double>& Atilde, const Eigen::VectorXd& C) {
    int N = grid.nx * grid.ny;
    
    // Compute diffusion term: Atilde * C
    Eigen::VectorXd rhs = Atilde * C;

    // Add chemical reaction terms if enabled
    if (grid.reactions) {
        for (int i = 0; i < N; ++i) {
            double c1 = C[i];
            double c2 = C[i + N];
            double c1_sq = c1 * c1;

            double omega1 = grid.Ca + c1_sq * c2 - grid.Cb * c1 - c1;
            double omega2 = grid.Cb * c1 - c1_sq * c2;

            rhs[i] += omega1;
            rhs[i + N] += omega2;
        }
    }

    return rhs;
}

Eigen::VectorXd computeResidual(const GridInfo& grid, const Eigen::SparseMatrix<double>& Atilde, const Eigen::VectorXd& C_new, const Eigen::VectorXd& C_old) {
    int N = grid.nx * grid.ny;
    
    // g = C_new - C_old - dt * (Atilde * C_new)
    Eigen::VectorXd g = C_new - C_old - grid.dt * (Atilde * C_new);
    
    // Add reaction term contribution if enabled
    if (grid.reactions) {
        for (int i = 0; i < N; ++i) {
            double c1 = C_new[i];
            double c2 = C_new[i + N];
            double c1_sq = c1 * c1;
            
            double omega1 = grid.Ca + c1_sq * c2 - grid.Cb * c1 - c1;
            double omega2 = grid.Cb * c1 - c1_sq * c2;
            
            g[i] -= grid.dt * omega1;
            g[i + N] -= grid.dt * omega2;
        }
    }
    
    return g;
}

Eigen::SparseMatrix<double> buildJacobian(const GridInfo& grid, const Eigen::SparseMatrix<double>& Atilde, const Eigen::VectorXd& C) {
    int N = grid.nx * grid.ny;
    Eigen::SparseMatrix<double> J(2 * N, 2 * N);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(Atilde.nonZeros() + 4 * N);
    
    // Add contribution of: I - dt * Atilde
    for (int k = 0; k < Atilde.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Atilde, k); it; ++it) {
            int r = it.row();
            int c = it.col();
            double val = -grid.dt * it.value();
            if (r == c) {
                val += 1.0;
            }
            triplets.push_back(Eigen::Triplet<double>(r, c, val));
        }
    }
    
    // Add local reaction derivative contributions if enabled
    if (grid.reactions) {
        for (int p = 0; p < N; ++p) {
            double c1 = C[p];
            double c2 = C[p + N];
            
            double d_w1_dc1 = 2.0 * c1 * c2 - grid.Cb - 1.0;
            double d_w1_dc2 = c1 * c1;
            double d_w2_dc1 = grid.Cb - 2.0 * c1 * c2;
            double d_w2_dc2 = -c1 * c1;
            
            triplets.push_back(Eigen::Triplet<double>(p, p, -grid.dt * d_w1_dc1));
            triplets.push_back(Eigen::Triplet<double>(p, p + N, -grid.dt * d_w1_dc2));
            triplets.push_back(Eigen::Triplet<double>(p + N, p, -grid.dt * d_w2_dc1));
            triplets.push_back(Eigen::Triplet<double>(p + N, p + N, -grid.dt * d_w2_dc2));
        }
    }
    
    J.setFromTriplets(triplets.begin(), triplets.end());
    return J;
}
