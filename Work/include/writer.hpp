/*
Write the current concentration fields C1 and C2 to a
modern XML VTK structured grid file (.vts)
*/

#ifndef WRITER_HPP
#define WRITER_HPP

#include "grid.hpp"
#include <Eigen/Core>
#include <string>

// Write the current concentration fields C1 and C2 to a
// modern XML VTK structured grid file (.vts)
void writeVTK(const std::string &filename, const GridInfo &grid,
              const Eigen::VectorXd &C);

#endif // WRITER_HPP
