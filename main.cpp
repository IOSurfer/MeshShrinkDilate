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
#include <vtkNIFTIImageReader.h>
#include <vtkNIFTIImageWriter.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkSmartPointer.h>

int main(int argc, char *argv[]) {
    if (argc < 9) {
        std::cerr << "Usage: " << argv[0]
                  << " InputStlFileName(.stl 0 for none)"
                  << " InputSdfFileName(.nii.gz 0 for none)"
                  << " OutputStlFileName(.stl)"
                  << " OutputSdfFileName(.nii.gz)"
                  << " DilateValue(mm)"
                  << " Spacing(mm)"
                  << " Padding(mm)"
                  << " SimplificationRate(0.0-1.0)"
                  << std::endl;
        return EXIT_FAILURE;
    }
    double dilate_value = atof(argv[5]);
    double spacing = atof(argv[6]);
    double padding = atof(argv[7]);
    double simplification_rate = atof(argv[8]);
    vtkSmartPointer<vtkImageData> image_data = vtkSmartPointer<vtkImageData>::New();

    if (strcmp(argv[1], "0") != 0) {
        vtkSmartPointer<vtkSTLReader> stl_reader = vtkSmartPointer<vtkSTLReader>::New();
        stl_reader->SetFileName(argv[1]);
        stl_reader->Update();

        vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
        decimate->SetInputData(stl_reader->GetOutput());
        decimate->SetTargetReduction(simplification_rate);
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
                    double test_point[3] = {(i + 0.5) * sdf.m_spacing_x + aabb.m_x_min, (j + 0.5) * sdf.m_spacing_y + aabb.m_y_min, (k + 0.5) * sdf.m_spacing_z + aabb.m_z_min};

                    double origin[3] = {(i + 0.5) * sdf.m_spacing_x + aabb.m_x_min, aabb.m_y_max, (k + 0.5) * sdf.m_spacing_z + aabb.m_z_min};
                    double closest_point[3];
                    double closest_point_dist2;
                    vtkIdType cell_id;
                    int sub_id;
                    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
                    cellLocator->FindClosestPoint(test_point, closest_point, cell, cell_id, sub_id, closest_point_dist2);
                    double distance = sqrt(closest_point_dist2);

                    // vtkSmartPointer<vtkIdList> cell_ids = vtkSmartPointer<vtkIdList>::New();
                    // cellLocator->FindCellsAlongLine(test_point, origin, 0.000001, cell_ids);

                    // if (distance <= spacing) {
                    //     distance = 0;
                    // }
                    // distance = (cell_ids->GetNumberOfIds() % 2 != 0) ? -distance : distance;
                    sdf.setValue(i, j, k, distance);
                }
            }
        }

        image_data->SetDimensions(range_x, range_y, range_z);
        image_data->SetSpacing(sdf.m_spacing_x, sdf.m_spacing_y, sdf.m_spacing_z);
        image_data->SetOrigin(aabb.m_x_min, aabb.m_y_min, aabb.m_z_min);
        image_data->AllocateScalars(VTK_DOUBLE, 1);

        for (int x = 0; x < range_x; ++x) {
            for (int y = 0; y < range_y; ++y) {
                for (int z = 0; z < range_z; ++z) {
                    int index = x * range_y * range_z + y * range_z + z;
                    double *pixel = static_cast<double *>(image_data->GetScalarPointer(x, y, z));
                    pixel[0] = sdf.m_volume_data[index];
                }
            }
        }
    } else if (strcmp(argv[2], "0") != 0) {
        vtkSmartPointer<vtkNIFTIImageReader> nifti_reader = vtkSmartPointer<vtkNIFTIImageReader>::New();
        nifti_reader->SetFileName(argv[2]);
        nifti_reader->Update();
        image_data = nifti_reader->GetOutput();
    } else {
        std::cout << "Invalid Arguments!"
                  << "\n";
        return EXIT_FAILURE;
    }

    vtkSmartPointer<vtkNIFTIImageWriter> nifti_writer = vtkSmartPointer<vtkNIFTIImageWriter>::New();
    nifti_writer->SetFileName(argv[4]);
    nifti_writer->SetInputData(image_data);
    nifti_writer->SetDescription("This is SDF");
    nifti_writer->SetQFac(1);
    nifti_writer->Write();

    vtkSmartPointer<vtkMarchingCubes> marchingCubes = vtkSmartPointer<vtkMarchingCubes>::New();
    marchingCubes->SetInputData(image_data);
    marchingCubes->SetValue(0, dilate_value);
    marchingCubes->Update();
    vtkSmartPointer<vtkPolyData> result_poly_data = marchingCubes->GetOutput();

    vtkSmartPointer<vtkSTLWriter> stl_writer = vtkSmartPointer<vtkSTLWriter>::New();
    stl_writer->SetFileName(argv[3]);
    stl_writer->SetInputData(result_poly_data);
    stl_writer->Write();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(result_poly_data);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(actor);
    renderer->SetBackground(0.1, 0.2, 0.4);

    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(800, 600);

    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}
