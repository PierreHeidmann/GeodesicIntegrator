#ifndef RENDER_HPP
#define RENDER_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Rk4vec.hpp"

template <typename data_t>
class render_black_hole
{
  public:
    void picture(int *red, int *green, int *blue, Vec3 center, double max_x,
                 double max_y, const int resolution,  const double alpha , const int start_ind,
                 int end_ind, const double TIME_MAXIMUM = 150.0,
                 const double DT = 0.1, const double T_START = 0);

    void render_circle(int *red, int *green, int *blue, double max_x,
                       double max_y,const int resolution,const int start_ind ,
                       int end_ind, double Radius = 2,
                       double thickness = 0.1);

    void render(int *red, int *green, int *blue, const int resolution,std::string file_name);
};

#endif
