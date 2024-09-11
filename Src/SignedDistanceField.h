#ifndef SIGNEDDISTANCEFIELD_H
#define SIGNEDDISTANCEFIELD_H

#include "AabbBox.h"
#include <string.h>

class SignedDistanceField {
  public:
    SignedDistanceField(AabbBox bounding_box, double spacing_x, double spacing_y, double spacing_z) : m_bouding_box(bounding_box), m_spacing_x(spacing_x), m_spacing_y(spacing_y), m_spacing_z(spacing_z) {
        m_range_x = (m_bouding_box.m_x_max - m_bouding_box.m_x_min) / m_spacing_x + 1;
        m_range_y = (m_bouding_box.m_y_max - m_bouding_box.m_y_min) / m_spacing_y + 1;
        m_range_z = (m_bouding_box.m_z_max - m_bouding_box.m_z_min) / m_spacing_z + 1;
        m_volume_data = new double[m_range_x * m_range_y * m_range_z];
    };

    ~SignedDistanceField() {
        if (m_volume_data != nullptr) {
            delete[] m_volume_data;
            m_volume_data = nullptr;
        }
    }

    void setValue(int x, int y, int z, double value) {
        m_volume_data[x * m_range_y * m_range_z + y * m_range_z + z] = value;
    }

    double getValue(int x, int y, int z) {
        return m_volume_data[x * m_range_y * m_range_z + y * m_range_z + z];
    }

    int getIndex(int x, int y, int z) {
        return x * m_range_y * m_range_z + y * m_range_z + z;
    }

    void getPosition(int x, int y, int z, double *points) {
        points[0] = (x + 0.5) * m_spacing_x + m_bouding_box.m_x_min;
        points[1] = (y + 0.5) * m_spacing_y + m_bouding_box.m_y_min;
        points[2] = (z + 0.5) * m_spacing_z + m_bouding_box.m_z_min;
    }

    AabbBox m_bouding_box;
    double m_spacing_x, m_spacing_y, m_spacing_z;
    int m_range_x, m_range_y, m_range_z;
    double *m_volume_data{nullptr};
};

#endif