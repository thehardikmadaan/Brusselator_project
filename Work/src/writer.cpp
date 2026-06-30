#include "writer.hpp"
#include <fstream>
#include <iostream>

void writeVTK(const std::string& filename, const GridInfo& grid, const Eigen::VectorXd& C) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing!\n";
        return;
    }

    int N = grid.nx * grid.ny;

    // Write legacy VTK header
    file << "# vtk DataFile Version 3.0\n";
    file << "Brusselator Concentration Output\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << grid.nx << " " << grid.ny << " 1\n";
    file << "ORIGIN 0.0 0.0 0.0\n";
    file << "SPACING " << grid.dx << " " << grid.dy << " 1.0\n";
    file << "POINT_DATA " << N << "\n";

    // Write scalar field for Species 1 (C1 - activator)
    file << "SCALARS C1 double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < N; ++i) {
        file << C[i] << "\n";
    }

    // Write scalar field for Species 2 (C2 - inhibitor)
    file << "SCALARS C2 double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < N; ++i) {
        file << C[i + N] << "\n";
    }

    file.close();
}
