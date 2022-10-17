#ifndef TOPOLOGICALSTAR5D_HPP
#define TOPOLOGICALSTAR5D_HPP

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "DimensionDefinitions.hpp"
#include "Rk4vec.hpp"
#include "TensorAlgebra.hpp"
#include "tensor.hpp"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

using namespace std;

// Do not change Dimension without talking to me :) you need to do more than just changing dimensions 
#define DIM 5

class Topological_Star5d
{

  public:
    const int dim = DIM;

    static tensor<2, double, DIM> get_metric(double rb, double rs, double x, double y, double z, double t,double theta = 0 /*temporary extra dimension that defaults to zero*/ );

    static tensor<3, double, DIM> get_metric_deriv(double rb, double rs, double x, double y, double z,
                                              double t,double theta = 0 /*temporary extra dimension that defaults to zero*/);

    static int eval_diff_eqn(double t, const double y[], double f[],
                             void *params);

    static gsl_matrix* invert_a_matrix(gsl_matrix *matrix);
    static tensor<2, double,DIM> invert_a_matrix( tensor<2, double,DIM> g);

};

#endif