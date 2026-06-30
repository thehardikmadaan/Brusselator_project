#ifndef WRITER_HPP
#define WRITER_HPP

#include "grid.hpp"
#include <Eigen/Core>
#include <string>

// Write the current concentration fields C1 and C2 to a legacy ASCII VTK structured points file
void writeVTK(const std::string& filename, const GridInfo& grid, const Eigen::VectorXd& C);

#endif // WRITER_HPP
