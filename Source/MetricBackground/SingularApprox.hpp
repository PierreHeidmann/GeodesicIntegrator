#ifndef SINGULARAPPROX_HPP
#define SINGULARAPPROX_HPP

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "DimensionDefinitions.hpp"
#include "Rk4vec.hpp"
#include "TensorAlgebra.hpp"
#include "tensor.hpp"

using namespace std;

class Singular_Approx
{

  public:
    static tensor<2, double> get_metric(double M, double x, double y, double z,double t);

    static tensor<3, double> get_metric_deriv(double M, double x, double y,
                                              double z, double t);


    static int eval_diff_eqn(double t, const double y[], double f[],
                             void *params);

};

#endif
