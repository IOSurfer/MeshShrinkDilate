#ifndef AABBBOX_H
#define AABBBOX_H

class AabbBox {
  public:
    AabbBox(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max) : m_x_min(x_min),
                                                                                                  m_x_max(x_max),
                                                                                                  m_y_min(y_min),
                                                                                                  m_y_max(y_max),
                                                                                                  m_z_min(z_min),
                                                                                                  m_z_max(z_max){

                                                                                                  };
    double m_x_min;
    double m_x_max;
    double m_y_min;
    double m_y_max;
    double m_z_min;
    double m_z_max;
};

#endif