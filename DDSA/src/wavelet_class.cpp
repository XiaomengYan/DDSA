#include "wavelet_class.hpp"

arma::mat wavelet_class::compute_wavelet_mat(arma::vec const &lambda){
  arma::mat res_wavelet_mat(lambda.n_elem, n_coeff), res_spline_mat(lambda.n_elem, n_spline_vec);
  arma::vec tmp_wavelet, tmp_spline;
  for( int i = 0; i < lambda.n_elem; i++){
    tmp_wavelet = compute_wavelet_vec(lambda(i));
    res_wavelet_mat.row(i) = tmp_wavelet.t();
    tmp_spline = compute_spline_vec(lambda(i));
    res_spline_mat.row(i) = tmp_spline.t();
  }
  return join_rows(res_spline_mat, res_wavelet_mat);
}

arma::vec wavelet_class::compute_wavelet_vec(double x){
  int k, j, tmp_n = 0, numK;
  arma::vec  coeff(n_coeff);
  for( k = 0 ; k < k_level; k++){
    numK = J_max(k) - J_min(k) + 1;
    for( j = 0; j < numK ; j++){
      coeff(j + tmp_n) = gen_mother_wavelet(x, k, j+J_min(k));
    }
    tmp_n += numK;
  }
  return coeff;
}

inline arma::vec wavelet_class::compute_spline_vec(double x){
  arma::vec res = spline_vec_obj.compute_spline(x);
  return res;
}


int wavelet_class::set_parameter(double lambda_min_, double lambda_max_,
                                 double deltalambda_, int k_level_){
  lambda_min = lambda_min_;
  lambda_max = lambda_max_;
  deltalambda = deltalambda_;
  k_level = k_level_;
  int k, J_maxC, J_minC;
  J_max = vec(k_level);
  J_min = vec(k_level);
  n_coeff = 0;
  for( k = 0 ; k < k_level; k++){
    J_maxC = ceil(pow(2,1.0*k)*lambda_max/deltalambda) - 1;
    J_minC = floor(pow(2,1.0*k)*lambda_min/deltalambda) - 6;
    J_max(k) = J_maxC;
    J_min(k) = J_minC;
    n_coeff += J_maxC - J_minC + 1;
  }
  n_knots = (lambda_max - lambda_min)/deltalambda - 1;
  n_spline_vec = spline_vec_obj.set_range(lambda_min, lambda_max, n_knots, 4);
  return n_coeff + n_spline_vec;
}




double wavelet_class::gen_mother_wavelet(double x, int k, int j){
  double y = pow(2, 1.0 * k)*x/deltalambda - j;
  return pow(2.0,1.0*k/2.0)*mother_wavelet(y);
}


double wavelet_class::mother_wavelet(double x){
  double coeff[7]={0.0001984127,0.0238095238, 0.2363095238,
                   0.4793650794, 0.2363095238, 0.0238095238,
                   0.0001984127};
  int j;
  double tmp, res,sign;
  res = 0;
  sign = 1;
  for(j=0;j<=6;j++){
    tmp = compute_splineD1(2*x-j);
    res += sign*coeff[j]*tmp;
    sign *= -1;
  }
  return res/8.0;
}

double wavelet_class::compute_splineD1(double x){
  if(x < 0 || x > 8)
    return 0.0;
  int m = 4, nbreaks = 23, n;
  n = 2*m  + nbreaks - 2;

  RcppGSL::matrix<double> coeff(n, m+1);
  gsl_bspline_deriv_eval(x, m, coeff, bw);

  return coeff(14, m);				// return vector
}



RCPP_MODULE(wavelet_class){
  Rcpp::class_<wavelet_class>("wavelet_class")
  .constructor()
  .method("compute_wavelet_mat", &wavelet_class::compute_wavelet_mat)
  .method("compute_wavelet_vec", &wavelet_class::compute_wavelet_vec)
  .method("compute_spline_vec", &wavelet_class::compute_spline_vec)
  .method("gen_mother_wavelet", &wavelet_class::gen_mother_wavelet)
  .method("set_parameter", &wavelet_class::set_parameter)
  .method("mother_wavelet", &wavelet_class::mother_wavelet);
}

