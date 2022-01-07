#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_bspline.h>
#include <cmath>
#include <vector>
#include "spline_class.hpp"


using namespace std;
using namespace arma;


class wavelet_class{

public:

  double compute_splineD1(double x);

  double mother_wavelet(double x);

  double gen_mother_wavelet(double x, int k, int j);

  arma::vec compute_wavelet_vec(double x);

  int set_parameter(double lambda_min_, double lambda_max_, double deltalambda_, int k_level_);

  arma::mat compute_wavelet_mat(arma::vec const &lambda);

  arma::vec compute_spline_vec(double x);

  wavelet_class(){
    int m = 4, nbreaks = 23, n;
    n = 2*m + nbreaks - 2;
    bw = gsl_bspline_alloc(2*m, nbreaks);
    gsl_bspline_knots_uniform(-7,15,bw);
  }

  ~wavelet_class(){
    gsl_bspline_free(bw);
  }

private:
  arma::vec lambda;
  gsl_bspline_workspace *bw;
  double lambda_min, lambda_max, deltalambda;
  int n_knots;
  int n_coeff, k_level, n_spline_vec;
  arma::vec J_max, J_min;
  spline_class spline_vec_obj;
};
