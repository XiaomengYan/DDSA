#include "spline_class.hpp"

//Function to evaluate the number of spline functions
int spline_class::set_range(double t_min_, double t_max_, int n_knots_, int m_order_){
  n_knots = n_knots_;
  t_min = t_min_;
  t_max = t_max_;
  m_order = m_order_;
  n_coeff = m_order + n_knots - 2;
  bw = gsl_bspline_alloc(m_order, n_knots);
  gsl_bspline_knots_uniform(t_min,t_max,bw);
  arma_coeff = arma::vec(n_coeff);
  return n_coeff;
}

arma::vec spline_class::compute_spline(double t){
  RcppGSL::vector<double> coeff(n_coeff);
  gsl_bspline_eval(t,coeff,bw);
  for(int i = 0 ; i < n_coeff; i++){
    arma_coeff(i) = coeff[i];
  }
  return arma_coeff;
}

RCPP_MODULE(spline_class){
  Rcpp::class_<spline_class>("spline_class")
  .constructor()
  .method("set_range", &spline_class::set_range)
  .method("compute_spline", &spline_class::compute_spline);
}
