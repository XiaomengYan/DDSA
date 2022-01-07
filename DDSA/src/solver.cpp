#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat GrLassoCpp(arma::mat Z, double lambdainvL){
  int p = Z.n_rows;
  arma::vec flagVec(2, fill::zeros);
  for(int i = 0; i < p; i++){
    flagVec(1) = 1.0 - lambdainvL/norm(Z.row(i), 2);
    Z.row(i) = flagVec.max()*Z.row(i);
  }
  return(Z);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat ProxGrLassoCpp(double lambda1, double lambda2, arma::mat WtY, List WtWList, List BtBList, arma::mat Omega){
  int K = WtWList.size();
  std::vector<arma::mat> WtW(K);
  std::vector<arma::mat> BtB(K);
  for(int k = 0; k < K; k++){
    SEXP Wl = WtWList[k];
    SEXP Bl = BtBList[k];
    Rcpp::NumericMatrix WRcpp(Wl);
    Rcpp::NumericMatrix BRcpp(Bl);
    WtW[k] = arma::mat(WRcpp.begin(), WRcpp.nrow(), WRcpp.ncol());
    BtB[k] = arma::mat(BRcpp.begin(), BRcpp.nrow(), BRcpp.ncol());
  }
  //std::cout << "Initialization Done!" << std::endl;
  int p = WtY.n_rows, q = WtY.n_cols;
  double diff = 1000.0;
  double ABSTOL = 0.01;
  int maxiter = 1000;
  int iter = 0;
  double L = 0.0;
  bool flag = true;

  // initialization
  arma::mat preX(p,q, fill::zeros); //x0
  arma::mat preY(preX), currX(preX), currY(preX); // y0,x1,y1
  double pret = 1.0, currt;
  double normCurrY = 0.0;

  while(flag && (iter < maxiter)){
    //std::cout << iter << std::endl;
    iter++;
    //L = pow(2000.0 + 1.0*iter, 0.8); //Kaepora Training
    //L = pow(500.0 + 1.0*iter, 0.8); // When the error is very small, this step size should
    L = pow(5000.0 + 1.0 * iter, 0.8);// SNEMO log training initial stepsize
    arma::mat tmpMat(p,q,fill::zeros);
    for(int k = 0; k < K; k ++){
      tmpMat += WtW[k]*preY*BtB[k];
    }
    currX = preY +  1.0/L * (1.0/K * (WtY - tmpMat) - lambda2 * preY * Omega);
    currX = GrLassoCpp(currX, lambda1/L);
    currt = (1.0 + sqrt(1.0 + 4.0*pow(pret,2.0)))/2.0;
    currY = currX + (pret - 1.0)/currt *(currX - preX);
    normCurrY = norm(currY,2);
    if(normCurrY == 0.0){
      iter = maxiter;
      std::cout << "All Elements in the Coefficient Matrix are ZEROS!" << std::endl;
    }else{
      diff = norm(currY - preY, "fro");
      if(diff > 10000.0){
        iter = maxiter;
        std::cout << "The Algorithm Does Not Converge!" << std::endl;
      }
      //Update
      if(diff < ABSTOL) flag = false;
      pret = currt;
      preY = currY;
      preX = currX;
    }
  }
  return(preY);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat MatGrLassoCpp(arma::mat Z, double etainL, int p, int q){
  int i;
  arma::vec flagVec(2, fill::zeros);
  for(i = 0; i < p ; i++){
    arma::mat subMat = Z.rows(span(i*q,(i+1)*q-1));
    double normMat = norm(subMat, "fro");
    if(normMat == 0.0){
      flagVec(1) = 0;
    }else{
      flagVec(1) = 1.0 - etainL/norm(subMat, "fro");
    }
    Z.rows(span(i*q,(i+1)*q-1)) = flagVec.max()*Z.rows(span(i*q,(i+1)*q-1));
  }
  return(Z);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double CV_PCTrainingErrCpp(arma::mat U, arma::vec TrainIndex, List SNeList, std::vector<arma::mat> AtA, std::vector<arma::mat> AtY, int p, int q){
  int S = TrainIndex.n_elem;
  int R = U.n_cols;
  int t = 0;
  arma::vec Vs(R,fill::zeros);
  arma::mat CoeffVec(U.n_rows, 1, fill::zeros);
  arma::vec Err(2000, fill::zeros); // create a large enough vector to stores the errors
  for(int s = 0; s < S;s++){
    int Index = TrainIndex[s];
    Environment SNeObj = SNeList[Index-1];
    Rcpp::List FluxList = SNeObj["FluxMean0List"];
    Rcpp::List WmatList = SNeObj["WmatList"];
    Rcpp::List BmatList = SNeObj["BmatList"];
    Vs = pinv(U.t()*AtA[Index-1]*U)*U.t()*AtY[Index-1];
    CoeffVec = U*Vs;
    arma::mat CoeffMat = arma::mat(CoeffVec.begin(), q, p).t();
    int N = SNeObj["N"];
    for(int n = 0; n<N;n++){
      arma::mat Flux = FluxList.at(n);
      arma::mat Wmat = WmatList.at(n);
      arma::mat Bmat = BmatList.at(n);
      Err[t] = pow(norm(Flux - Wmat*CoeffMat*Bmat, "fro"),2)/Flux.n_elem;
      t++;
    }

  }
  return(sum(Err)/t);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List ADAMCpp(arma::vec TrainIndex, int R, double eta1, double eta2, List SNeList, List AtAList, List AtYList, arma::mat Omega){
  int K = AtAList.size();
  std::vector<arma::mat> AtA(K);
  std::vector<arma::mat> AtY(K);

  for(int k = 0; k < K; k++){
    SEXP AtAl = AtAList[k];
    SEXP AtYl = AtYList[k];
    Rcpp::NumericMatrix AtARcpp(AtAl);
    Rcpp::NumericMatrix AtYRcpp(AtYl);
    AtA[k] = arma::mat(AtARcpp.begin(), AtARcpp.nrow(), AtARcpp.ncol());
    AtY[k] = arma::mat(AtYRcpp.begin(), AtYRcpp.nrow(), AtYRcpp.ncol());
  }
  //std::cout << "Initialization Done!" << std::endl;
  int q = Omega.n_cols;
  int p = AtA[0].n_rows/q;
  arma::vec Ivec(p,fill::ones);
  arma::mat IOmega = eta2*kron(diagmat(Ivec), Omega);
  double beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8;
  int iter = 0, n = 0;
  double diff = 1000.0;
  bool flag = true;
  double ABSTOL = 1e-6;
  int maxiter = 20000;
  arma::vec MseVec(200, fill::zeros);

  // Initialization
  arma::mat preU(p*q,R,fill::randn);
  arma::mat tmpQ, tmpR;
  qr_econ(tmpQ, tmpR, preU);
  preU = tmpQ;
  arma::mat preM(p*q, R, fill::zeros), preV(preM), currM(preM), currV(preM), currU(preM), currM_corr(currM), currV_corr(currV);

  //
  arma::mat Vs(R,1, fill::zeros);
  arma::mat Grad(p*q, R, fill::zeros);
  arma::mat tmpMat(R,R,fill::zeros);


  while(flag && iter<maxiter){
    iter++;
    //std::cout << iter << std::endl;
    arma::uvec index = randperm(TrainIndex.n_elem,1);
    int s_train = TrainIndex(index.at(0));
    double L = pow((100.0 + 1.0*iter), 0.98);
    Vs = pinv(preU.t()*AtA[s_train-1]*preU,1e-4)*preU.t()*AtY[s_train-1];
    //std::cout << Vs << std::endl;
    Grad = -1.0*(AtY[s_train-1] - AtA[s_train-1] * preU * Vs) * Vs.t() + IOmega * preU;
    //std::cout << norm(Grad,"fro") << std::endl;
    currM = beta1*preM + (1.0 - beta1) * Grad;
    //std::cout << norm(currM,"fro") << std::endl;
    currV = beta2 * preV + (1.0 - beta2) * pow(Grad,2);
    //std::cout << norm(currV,"fro") << std::endl;
    currM_corr = currM/(1.0 - pow(beta1, iter));
    //std::cout << norm(currM,"fro") << std::endl;
    currV_corr = currV/(1.0 - pow(beta2, iter));
    Grad = currM_corr/(sqrt(currV_corr) +  epsilon);
    tmpMat = preU.t() * Grad;
    Grad = Grad - 0.5 * preU * (tmpMat + tmpMat.t());
    currU = MatGrLassoCpp(preU - 1.0/L*Grad, eta1/L, p, q);
    qr_econ(tmpQ, tmpR, currU);
    currU = tmpQ;
    //std::cout << currU << std::endl;
    // double normU = norm(currU-preU, "fro");
    // std::cout << normU << std::endl;

    if(diff > 10000.0){
      iter = maxiter;
      std::cout << "The Algorithm Does Not Converge!" << std::endl;
    }
    // //Update
    if(iter%100 == 0){
      MseVec(n) = CV_PCTrainingErrCpp(currU, TrainIndex, SNeList, AtA, AtY, p, q);
      //std::cout << MseVec(n) << std::endl;
      if(n > 0) diff = abs(MseVec(n-1) - MseVec(n));
      if(diff< ABSTOL){flag = false;}
      n++;
    }
    //Update
    preU = currU;
    preM = currM;
    preV = currV;

  }
  arma::mat V(R, TrainIndex.n_elem, fill::zeros);
  for(int i = 0; i < TrainIndex.n_elem; i ++){
    int index = TrainIndex(i);
    V.col(i) = pinv(preU.t()*AtA[index-1]*preU)*preU.t()*AtY[index-1];
  }
  return Rcpp::List::create(Rcpp::Named("U0") = preU,Rcpp::Named("V0") = V,Rcpp::Named("MSE_train") = MseVec);
}






