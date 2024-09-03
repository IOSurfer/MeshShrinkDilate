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
    // Ensure a filename is provided
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " Filename(.stl)"
                  << " Dilate value"
                  << " Spacing(cm)"
                  << " Padding(cm)" << std::endl;
        return EXIT_FAILURE;
    }
    double dilate_value = atof(argv[2]);
    double spacing = atof(argv[3]);
    double padding = atof(argv[4]);

    // Create a reader for the STL file
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();

    // 使用 vtkDecimatePro 进行三角面简化
    vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
    decimate->SetInputData(reader->GetOutput());
    decimate->SetTargetReduction(1.0); // 将三角形数目减少50%
    decimate->PreserveTopologyOn();    // 保持拓扑结构
    decimate->Update();
    vtkSmartPointer<vtkPolyData> poly_data = decimate->GetOutput();

    if (!poly_data || poly_data->GetNumberOfPoints() == 0) {
        std::cerr << "Error: No data was read from the STL file." << std::endl;
        return EXIT_FAILURE;
    }

    vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
    cellLocator->SetDataSet(poly_data);
    cellLocator->BuildLocator();
    // vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    // implicitPolyDataDistance->SetInput(poly_data);

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
                // Example point to test
                std::array<double, 3> testPoint = {(i + 0.5) * sdf.m_spacing_x + aabb.m_x_min, (j + 0.5) * sdf.m_spacing_y + aabb.m_y_min, (k + 0.5) * sdf.m_spacing_z + aabb.m_z_min}; // Replace with your test point

                // Find the closest point on the mesh
                double closestPoint[3];
                double closestPointDist2;
                vtkIdType cellId;
                int subId;
                vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();

                cellLocator->FindClosestPoint(testPoint.data(), closestPoint, cell, cellId, subId, closestPointDist2);

                // Calculate the distance (already squared)
                double distance = sqrt(closestPointDist2);
                // double distance = implicitPolyDataDistance->EvaluateFunction(testPoint.data());
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

    // 3. 使用 Marching Cubes 算法提取等值面
    vtkSmartPointer<vtkMarchingCubes> marchingCubes = vtkSmartPointer<vtkMarchingCubes>::New();
    marchingCubes->SetInputData(imageData);
    marchingCubes->SetValue(0, dilate_value); // 设置等值面值

    // 4. 获取生成的多边形数据
    marchingCubes->Update();
    vtkSmartPointer<vtkPolyData> polyData = marchingCubes->GetOutput();

    // Create a mapper
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polyData);

    // Create an actor
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    // Create a renderer
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(actor);
    renderer->SetBackground(0.1, 0.2, 0.4); // Set background color to dark blue

    // Create a render window
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(800, 600);

    // Create a render window interactor
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Start the interaction
    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}
