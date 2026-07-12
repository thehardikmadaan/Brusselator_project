#include "grid.hpp"
#include "brusselator.hpp"
#include "integrator.hpp"
#include "writer.hpp"
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iomanip>
#include <Eigen/Core>

int main(int argc, char *argv[]) {
    // Determine which CLI format is being used           \\
    // 6-positional format: nx ny duration write_interval \\
    // Ca Cb [method] [reactions] (uses default dt=0.01)  \\
    // 7-positional format: nx ny duration dt             \\
    // write_interval Ca Cb [method] [reactions]          \\
    
    bool is_six_positional = false;
    
    if (argc == 7) {
        is_six_positional = true;
    } else if (argc == 8) {
        std::string arg7 = argv[7];
        if (arg7 == "rk4" || arg7 == "implicit" || arg7 == "explicit") {
            is_six_positional = true;
        }
    } else if (argc == 9) {
        std::string arg7 = argv[7];
        if (arg7 == "rk4" || arg7 == "implicit" || arg7 == "explicit") {
            is_six_positional = true;
        }
    } else if (argc < 7 || argc > 10) {
        std::cerr << "Usage (7 positional arguments):\n"
                  << "  " << argv[0] << " <nx> <ny> <duration> <dt> <write_interval> <Ca> <Cb> [method] [reactions]\n"
                  << "Usage (6 positional arguments, uses default dt=0.01):\n"
                  << "  " << argv[0] << " <nx> <ny> <duration> <write_interval> <Ca> <Cb> [method] [reactions]\n";
        return 1;
    }

    GridInfo grid;

    if (is_six_positional) {
        grid.nx = std::atoi(argv[1]);
        grid.ny = std::atoi(argv[2]);
        grid.duration = std::atof(argv[3]);
        grid.dt = 0.01; // default
        grid.write_interval = std::atof(argv[4]);
        grid.Ca = std::atof(argv[5]);
        grid.Cb = std::atof(argv[6]);
        
        if (argc >= 8) {
            grid.method = argv[7];
        }
        if (argc >= 9) {
            std::string react_str = argv[8];
            grid.reactions = (react_str == "on");
        }
    } else {
        grid.nx = std::atoi(argv[1]);
        grid.ny = std::atoi(argv[2]);
        grid.duration = std::atof(argv[3]);
        grid.dt = std::atof(argv[4]);
        grid.write_interval = std::atof(argv[5]);
        grid.Ca = std::atof(argv[6]);
        grid.Cb = std::atof(argv[7]);
        
        if (argc >= 9) {
            grid.method = argv[8];
        }
        if (argc >= 10) {
            std::string react_str = argv[9];
            grid.reactions = (react_str == "on");
        }
    }

    // Set other grid parameters                          \\
    grid.dx = 0.005;
    grid.dy = 0.005;
    grid.D1 = 1e-5;
    grid.D2 = 1e-6;

    std::cout << "Starting 2D Brusselator Simulation...\n";
    std::cout << "Grid: " << grid.nx << "x" << grid.ny << "\n";
    std::cout << "Domain: " << (grid.nx - 1) * grid.dx << "m x " << (grid.ny - 1) * grid.dy << "m\n";
    std::cout << "Method: " << (grid.method == "implicit" ? "Implicit Euler (Newton)" : (grid.method == "explicit" ? "Explicit Euler" : "Explicit RK4")) << "\n";
    std::cout << "Reactions: " << (grid.reactions ? "Enabled" : "Disabled") << "\n";
    std::cout << "Time step: dt = " << grid.dt << " s\n";

    // Initialize state                                   \\
    int num_points = grid.nx * grid.ny;
    Eigen::VectorXd C1 = Eigen::VectorXd::Constant(num_points, grid.Ca);
    Eigen::VectorXd C2 = Eigen::VectorXd::Constant(num_points, grid.Cb / grid.Ca);

    // Apply perturbations                                \\
    if (grid.nx > 30 && grid.ny > 30) {
        int p1_idx = 10 * grid.nx + 10;
        int p2_idx = 30 * grid.nx + 30;

        C1[p1_idx] *= 1.2;
        C2[p1_idx] *= 1.2;

        C1[p2_idx] *= 0.9;
        C2[p2_idx] *= 0.9;
        std::cout << "Applied initial perturbations.\n";
    }

    // Combine into monolithic state vector C             \\
    Eigen::VectorXd C(2 * num_points);
    C.head(num_points) = C1;
    C.tail(num_points) = C2;

    // Build discretization matrices                      \\
    std::cout << "Assembling sparse matrices...\n";
    Eigen::SparseMatrix<double> A = buildMatrixA(grid);
    Eigen::SparseMatrix<double> Atilde = buildMatrixAtilde(grid, A);

    // Create output directory                            \\
    #ifdef _WIN32
        system("mkdir output > nul 2>&1");
    #else
        system("mkdir -p output");
    #endif

    // Write initial state                                \\
    std::string init_filename = "output/output_0000.vts";
    writeVTK(init_filename, grid, C);

    double t = 0.0;
    int step = 0;
    int total_newton_its = 0;
    double last_write_time = 0.0;
    int file_counter = 1;

    std::cout << "Starting time integration...\n";
    while (t < grid.duration) {
        double current_dt = grid.dt;
        if (t + current_dt > grid.duration) {
            current_dt = grid.duration - t;
        }

        // We temporarily update grid.dt to match the     \\
        // actual step size taken                         \\
        GridInfo current_grid = grid;
        current_grid.dt = current_dt;

        if (grid.method == "implicit") {
            int prev_its = total_newton_its;
            stepImplicitEuler(current_grid, Atilde, C, total_newton_its);
            int step_its = total_newton_its - prev_its;

            t += current_dt;
            step++;

            if (step % 10 == 0 || t >= grid.duration) {
                std::cout << "Time: " << std::fixed << std::setprecision(3) << t << " s"
                          << " | Step: " << step 
                          << " | Newton Iterations: " << step_its << "\n";
            }
        } else if (grid.method == "explicit") {
            stepExplicitEuler(current_grid, Atilde, C);
            t += current_dt;
            step++;

            if (step % 100 == 0 || t >= grid.duration) {
                std::cout << "Time: " << std::fixed << std::setprecision(3) << t << " s"
                          << " | Step: " << step << "\n";
            }
        } else {
            stepRK4(current_grid, Atilde, C);
            t += current_dt;
            step++;

            if (step % 100 == 0 || t >= grid.duration) {
                std::cout << "Time: " << std::fixed << std::setprecision(3) << t << " s"
                          << " | Step: " << step << "\n";
            }
        }

        // Write output if interval reached               \\
        if (t - last_write_time >= grid.write_interval - 1e-9 || t >= grid.duration) {
            std::stringstream ss;
            ss << "output/output_" << std::setfill('0') << std::setw(4) << file_counter << ".vts";
            writeVTK(ss.str(), grid, C);
            last_write_time = t;
            file_counter++;
        }
    }

    std::cout << "Simulation completed successfully.\n";
    if (grid.method == "implicit") {
        std::cout << "Total Newton Iterations: " << total_newton_its << "\n";
    }
    return 0;
}
