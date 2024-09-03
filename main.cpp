#include "AabbBox.h"
#include "SignedDistanceField.h"
#include "string.h"
#include <vtkActor.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkDecimatePro.h>
#include <vtkGenericCell.h>
#include <vtkImageData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkMarchingCubes.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " Filename(.stl)"
                  << " Dilate value"
                  << " Spacing(mm)"
                  << " Padding(mm)" << std::endl;
        return EXIT_FAILURE;
    }
    double dilate_value = atof(argv[2]);
    double spacing = atof(argv[3]);
    double padding = atof(argv[4]);

    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();

    vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
    decimate->SetInputData(reader->GetOutput());
    decimate->SetTargetReduction(1.0);
    decimate->PreserveTopologyOn();
    decimate->Update();
    vtkSmartPointer<vtkPolyData> poly_data = decimate->GetOutput();

    if (!poly_data || poly_data->GetNumberOfPoints() == 0) {
        std::cerr << "Error: No data was read from the STL file." << std::endl;
        return EXIT_FAILURE;
    }

    vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
    cellLocator->SetDataSet(poly_data);
    cellLocator->BuildLocator();

    double bounding_box[6];
    poly_data->GetBounds(bounding_box);
    AabbBox aabb(bounding_box[0] - 5, bounding_box[1] + 5, bounding_box[2] - 5, bounding_box[3] + 5, bounding_box[4] - 5, bounding_box[5] + 5);
    SignedDistanceField sdf(aabb, spacing, spacing, spacing);
    std::cout << aabb.m_x_min << " " << aabb.m_x_max << " " << aabb.m_y_min << " " << aabb.m_y_max << " " << aabb.m_z_min << " " << aabb.m_z_max << "\n";
    int range_x, range_y, range_z;
    range_x = sdf.m_range_x;
    range_y = sdf.m_range_y;
    range_z = sdf.m_range_z;
    for (int i = 0; i < range_x; i++) {
        for (int j = 0; j < range_y; j++) {
            for (int k = 0; k < range_z; k++) {
                std::array<double, 3> testPoint = {(i + 0.5) * sdf.m_spacing_x + aabb.m_x_min, (j + 0.5) * sdf.m_spacing_y + aabb.m_y_min, (k + 0.5) * sdf.m_spacing_z + aabb.m_z_min}; // Replace with your test point

                double closestPoint[3];
                double closestPointDist2;
                vtkIdType cellId;
                int subId;
                vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();

                cellLocator->FindClosestPoint(testPoint.data(), closestPoint, cell, cellId, subId, closestPointDist2);

                double distance = sqrt(closestPointDist2);
                sdf.setValue(i, j, k, distance);
            }
        }
    }

    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(range_x, range_y, range_z);
    imageData->AllocateScalars(VTK_DOUBLE, 1);

    for (int x = 0; x < range_x; ++x) {
        for (int y = 0; y < range_y; ++y) {
            for (int z = 0; z < range_z; ++z) {
                int index = x * range_y * range_z + y * range_z + z;
                double *pixel = static_cast<double *>(imageData->GetScalarPointer(x, y, z));
                pixel[0] = sdf.m_volume_data[index];
            }
        }
    }

    vtkSmartPointer<vtkMarchingCubes> marchingCubes = vtkSmartPointer<vtkMarchingCubes>::New();
    marchingCubes->SetInputData(imageData);
    marchingCubes->SetValue(0, dilate_value);

    marchingCubes->Update();
    vtkSmartPointer<vtkPolyData> polyData = marchingCubes->GetOutput();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polyData);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(actor);
    renderer->SetBackground(0.1, 0.2, 0.4); // Set background color to dark blue

    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(800, 600);

    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}
