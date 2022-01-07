#include <RcppArmadillo.h>
#include <gsl/gsl_bspline.h>
#include <RcppGSL.h>
#include <cmath>
#include <vector>

using namespace std;
using namespace arma;

class spline_class{
public:

  int set_range(double t_min_, double t_max_, int n_kots_, int m_order_);

  arma::vec compute_spline(double t); //spline evaluation

  ~spline_class(){
    gsl_bspline_free(bw);
  }

  private:
    arma::vec arma_coeff;
    gsl_bspline_workspace *bw;
    int n_coeff, n_knots, m_order;
    double t_min, t_max;
    };
