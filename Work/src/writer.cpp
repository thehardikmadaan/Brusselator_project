#include "writer.hpp"
#include <iostream>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>

void writeVTK(const std::string& filename, const GridInfo& grid, const Eigen::VectorXd& C) {
    int N = grid.nx * grid.ny;

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            points->InsertNextPoint(i * grid.dx, j * grid.dy, 0.0);
        }
    }

    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
    structuredGrid->SetDimensions(grid.nx, grid.ny, 1);
    structuredGrid->SetPoints(points);

    vtkSmartPointer<vtkFloatArray> c1Array = vtkSmartPointer<vtkFloatArray>::New();
    c1Array->SetName("C1");
    c1Array->SetNumberOfComponents(1);
    c1Array->SetNumberOfTuples(N);

    vtkSmartPointer<vtkFloatArray> c2Array = vtkSmartPointer<vtkFloatArray>::New();
    c2Array->SetName("C2");
    c2Array->SetNumberOfComponents(1);
    c2Array->SetNumberOfTuples(N);

    for (int i = 0; i < N; ++i) {
        c1Array->SetValue(i, C[i]);
        c2Array->SetValue(i, C[i + N]);
    }

    structuredGrid->GetPointData()->AddArray(c1Array);
    structuredGrid->GetPointData()->AddArray(c2Array);

    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(structuredGrid);
    writer->SetCompressorTypeToZLib();
    writer->Write();
}
