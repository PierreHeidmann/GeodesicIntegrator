#ifndef OSCILLOTON_HPP
#define OSCILLOTON_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "DimensionDefinitions.hpp"
#include "Rk4vec.hpp"
#include "TensorAlgebra.hpp"
#include "tensor.hpp"

#define column_max 8
#define row_max 2600
#define spacing 0.01

class Oscilloton
{
  public:
    static double a202[row_max][8]; // a_j
    static double c202[row_max][8]; // c_j
    static double m_omega;

    Oscilloton();

    static int eval_diff_eqn(double t, const double y[], double f[],
                             void *params);

    static double calculate_norm(Vec3 v, double M = 1);

    // Change the vt component to set norm to any value
    static Vec3 set_norm(Vec3 v, double norm_val = 0, double M = 1);

    static double get_a202(double rr, int component);

    static double get_c202(double rr, int component);

  private:
    static tensor<2, double> get_metric(double M, double x, double y, double z,
                                        double t);

    static tensor<3, double> get_metric_deriv(double M, double x, double y,
                                              double z, double t);

    static tensor<3, double> get_chris(tensor<2, double> g_UU,
                                       tensor<3, double> dg);
};

#endif
