##############################################
## Step1 : Data Load In
##############################################

LoadData = function(DataDir, PhaseRange, WavelengthRange){
  PhaseMin = min(PhaseRange)
  PhaseMax = max(PhaseRange)
  WavelengthMin = min(WavelengthRange)
  WavelengthMax = max(WavelengthRange)
  WaveletObj = R_wavelet_class$new(WavelengthMin,WavelengthMax)
  SplineObj = R_spline_class$new(PhaseMin,PhaseMax)
  Omega = SplineObj$Omega
  SNeFiles = list.files(DataDir)
  S = length(SNeFiles)
  SNeList = list()
  for (s in 1:S) {
      Name = paste0("SN",unlist(strsplit(SNeFiles[s],split = ".dat")))
      SNedf = read.table(paste0(DataDir,SNeFiles[s]),header = TRUE)
      SNeList[[s]] = sn_class$new(SNedf = SNedf,SNe = Name,WaveletObj = WaveletObj,SplineObj = SplineObj)
    }
  return(list("SNeList" =SNeList,"Omega" = Omega))
}


###############################################
## Step2: Mean Training
###############################################


####################################################
# Mean Training Useful functions
#####################################################
TrainMeanPrepare = function(TrainIndex,SNeList){
  p0 = dim(SNeList[[1]]$WmatList[[1]])[2]
  q0 = dim(SNeList[[1]]$BmatList[[1]])[1]
  WtY = matrix(0, nrow = p0, ncol = q0)
  WtWList = list()
  BtBList = list()
  K = 0
  for (s_train in TrainIndex) {
    SNeObj = SNeList[[s_train]]
    N = SNeObj$N
    for (n in 1:N) {
      K = K+1
      tmpWmat = SNeObj$WmatList[[n]]/SNeObj$ErrList[[n]]
      WtY = WtY + crossprod(tmpWmat, SNeObj$FluxList[[n]]/SNeObj$ErrList[[n]])%*%t(SNeObj$BmatList[[n]])/length(SNeObj$FluxList[[n]])
      WtWList[[K]] = crossprod(tmpWmat)/length(SNeObj$FluxList[[n]])
      BtBList[[K]] = tcrossprod(SNeObj$BmatList[[n]])
    }
  }
  return(list("WtY" = WtY, "WtWList" = WtWList,"BtBList" = BtBList))
}

CV.MeanTrainingErr = function(Theta0, TestIndex, SNeList){
  Err = NULL
  for (s_test in TestIndex) {
    SNeObj = SNeList[[s_test]]
    N = SNeObj$N
    for (n in 1:N) {
      Err = c(Err, norm(SNeObj$FluxList[[n]] -  SNeObj$WmatList[[n]]%*%Theta0%*%SNeObj$BmatList[[n]],"2")^2/length(SNeObj$FluxList[[n]]))
    }
  }
  return(mean(Err))
}


CV.MeanTraining = function(lambda1Seq, lambda2Seq, SNeList, Omega, KFold = 5,ifParallel, NoCore){
  S = length(SNeList)
  CVgroup = createFolds(1:S, k = KFold)
  CVMat = matrix(0,nrow = length(lambda1Seq),ncol = length(lambda2Seq))
  if(ifParallel){
    cl <- makeCluster(NoCore)
    registerDoParallel(cl)
    ErrMatList <- foreach (k = 1:KFold) %dopar%  {
      #Rcpp::sourceCpp("../src/CoreFuns.cpp")
      cat(paste0("Fold = ", k))
      TestIndex = CVgroup[[k]]
      TrainIndex = setdiff(1:S, TestIndex)
      TrainPreList = TrainMeanPrepare(TrainIndex, SNeList)
      WtY = TrainPreList$WtY
      WtWList = TrainPreList$WtWList
      BtBList = TrainPreList$BtBList
      ErrMat = matrix(0,nrow = length(lambda1Seq),ncol = length(lambda2Seq))
      for (i in 1:length(lambda1Seq)) {
        lambda1 = lambda1Seq[i]
        for (j in 1:length(lambda2Seq)) {
          lambda2 = lambda2Seq[j]
          Theta0 = ProxGrLassoCpp(lambda1 = lambda1,lambda2 = lambda2,WtY = WtY,WtWList = WtWList,BtBList = BtBList,Omega = Omega)
          ErrMat[i,j] = CV.MeanTrainingErr(Theta0 = Theta0,TestIndex = TestIndex,SNeList = SNeList)
        }
      }
      return(ErrMat)
    }
    stopCluster(cl)
    for (k in 1:KFold) {
      CVMat = CVMat + ErrMatList[[k]]
    }
  }
  else{
    for (k in 1:KFold) {
      cat(paste0("Fold = ", k))
      TestIndex = CVgroup[[k]]
      TrainIndex = setdiff(1:S, TestIndex)
      TrainPreList = TrainMeanPrepare(TrainIndex, SNeList)
      WtY = TrainPreList$WtY
      WtWList = TrainPreList$WtWList
      BtBList = TrainPreList$BtBList
      ErrMat = matrix(0,nrow = length(lambda1Seq),ncol = length(lambda2Seq))
      for (i in 1:length(lambda1Seq)) {
        lambda1 = lambda1Seq[i]
        for (j in 1:length(lambda2Seq)) {
          lambda2 = lambda2Seq[j]
          Theta0 = ProxGrLassoCpp(lambda1,lambda2,Omega = Omega,WtY = WtY,WtWList = WtWList,BtBList = BtBList)
          ErrMat[i,j] = CV.MeanTrainingErr(Theta0 = Theta0,TestIndex = TestIndex,SNeList = SNeList)
        }
      }
      CVMat = CVMat + ErrMat
    }
  }
  BestPara = which(CVMat == min(CVMat), arr.ind = TRUE)
  lambda1 = lambda1Seq[BestPara[1]]
  lambda2 = lambda2Seq[BestPara[2]]
  return(list("CVScoreMat" = CVMat/KFold, "lambda1" = lambda1, "lambda2" = lambda2))
}


MeanTraining = function(lambda1, lambda2, SNeList, Omega){
  S = length(SNeList)
  TrainPreList = TrainMeanPrepare(1:S, SNeList)
  WtY = TrainPreList$WtY
  WtWList = TrainPreList$WtWList
  BtBList = TrainPreList$BtBList
  Theta0 = ProxGrLassoCpp(lambda1,lambda2,Omega = Omega,WtY = WtY,WtWList = WtWList,BtBList = BtBList) ## Use RCpp Code
  return(Theta0)
}

########################################################
## Step3: Principal Component Training
#########################################################

#######################################
## Train Principal Component Userful Functions
#################################
TrainPCPrepare = function(Theta0, SNeList){
  p0 = dim(Theta0)[1]
  q0 = dim(Theta0)[2]
  # Remove the mean from the training set
  S = length(SNeList)
  for (s in 1:S) {
    SNeObj = SNeList[[s]]
    SNeObj$RemoveMean(Theta0 = Theta0)
  }# create a new SNeList, don't forget to output
  ## Step 2: Calculate global list
  # We use A to represent W\otimes B.t()
  AtAList = list()
  AtYList = list()
  for (s in 1:S) {
    SNeObj = SNeList[[s]]
    N = SNeObj$N
    AtAmat = matrix(0,p0*q0,p0*q0)
    AtYmat = matrix(0,p0*q0,1)
    for (n in 1:N) {
      tmpWmat = SNeObj$WmatList[[n]]/SNeObj$ErrList[[n]]
      #Amat = kronecker(tmpWmat,t(SNeObj$BmatList[[n]]))
      AtAmat = AtAmat + kronecker(crossprod(tmpWmat), tcrossprod(SNeObj$BmatList[[n]]))/length(SNeObj$FluxList[[n]])
      AtYmat = AtYmat + kronecker(t(tmpWmat), SNeObj$BmatList[[n]]) %*%(SNeObj$FluxMean0List[[n]]/SNeObj$ErrList[[n]])/length(SNeObj$FluxList[[n]])
    }
    AtAList[[s]] = AtAmat/N
    AtYList[[s]] = AtYmat/N
  }
  return(list("SNeList" = SNeList, "AtAList" = AtAList, "AtYList" = AtYList))
}


CV.PCTrainingErr = function(U, TestIndex, SNeList, AtAList, AtYList){
  Err = NULL
  p0 = dim(SNeList[[1]]$WmatList[[1]])[2]
  q0 = dim(SNeList[[1]]$BmatList[[1]])[1]
  for (s_test in TestIndex) {
    SNeObj = SNeList[[s_test]]
    N = SNeObj$N
    Vs = ginv(t(U)%*%AtAList[[s_test]]%*%U)%*%t(U)%*%AtYList[[s_test]]
    CoeffMat = matrix(rowSums(sweep(U,MARGIN = 2,STATS = Vs, FUN = "*")), nrow = p0, ncol = q0, byrow = TRUE)
    for (n in 1:N) {
      Err = c(Err, norm(SNeObj$FluxMean0List[[n]]-SNeObj$WmatList[[n]]%*%CoeffMat%*%SNeObj$BmatList[[n]],"2")^2/length(SNeObj$FluxMean0List[[n]]))
    }
  }
  return(mean(Err))
}



CV.PCTraining = function(Theta0, RSeq, eta1Seq, eta2Seq,SNeList,Omega,KFold = 5,ifParallel = TRUE, NoCore = 5){
  TrainPCPreList = TrainPCPrepare(Theta0 = Theta0,SNeList = SNeList)
  SNeList = TrainPCPreList$SNeList
  AtAList = TrainPCPreList$AtAList
  AtYList = TrainPCPreList$AtYList
  S = length(SNeList)

  CVgroup = createFolds(1:S, k = KFold)
  CVErrList = list()
  if(ifParallel){
    cl <- makeCluster(NoCore)
    registerDoParallel(cl)
    CVErrList <- foreach(r = 1:length(RSeq),.packages = c("FpcaSED"),.noexport = c("ADAMCpp", "TrainPCPrepare","CV.PCTrainingErr")) %dopar%{
      Rcpp::sourceCpp("../FpcaSED/src/CoreFuns.cpp")
      cat(paste0("R = ", r,"\n"))
      R = RSeq[r]
      CVMat = matrix(0, nrow = length(eta1Seq), ncol = length(eta2Seq))
      for (k in 1:KFold) {
        TestIndex = CVgroup[[k]]
        TrainIndex = setdiff(1:S, TestIndex)
        ErrMat = matrix(0,nrow = length(eta1Seq),ncol = length(eta2Seq))
        for (i in 1:length(eta1Seq)) {
          eta1 = eta1Seq[i]
          for (j in 1:length(eta2Seq)) {
            eta2 = eta2Seq[j]
            ADAMRes = ADAMCpp(TrainIndex,R,eta1, eta2, SNeList = SNeList,AtAList = AtAList,AtYList = AtYList,Omega = Omega)
            ErrMat[i,j] = CV.PCTrainingErr(U = ADAMRes$U0,TestIndex = TestIndex,SNeList = SNeList,AtAList = AtAList,AtYList = AtYList)
          } # eta2
        } # eta1
        CVMat = CVMat + ErrMat
      }
      return(list("R" = R,"CVMat" = CVMat/KFold))
    }# loop
    stopCluster(cl)
  }else{
    for (r in 1:length(RSeq)) {
      cat(paste0("R = ", RSeq[r],"\n"))
      R = RSeq[r]
      CVMat = matrix(0, nrow = length(eta1Seq), ncol = length(eta2Seq))
      for (k in 1:KFold) {
        TestIndex = CVgroup[[k]]
        TrainIndex = setdiff(1:S, TestIndex)
        ErrMat = matrix(0,nrow = length(eta1Seq),ncol = length(eta2Seq))
        for (i in 1:length(eta1Seq)) {
          eta1 = eta1Seq[i]
          for (j in 1:length(eta2Seq)) {
            eta2 = eta2Seq[j]
            ADAMRes = ADAMCpp(TrainIndex,R,eta1, eta2, SNeList = SNeList,AtAList = AtAList,AtYList = AtYList,Omega = Omega)
            ErrMat[i,j] = CV.PCTrainingErr(U = ADAMRes$U0,TestIndex = TestIndex,SNeList = SNeList,AtAList = AtAList,AtYList = AtYList)
          } # eta2
        } # eta1
        CVMat = CVMat + ErrMat
      }
      CVErrList[[r]] = list("R" = R, "CVMat" = CVMat/KFold)
    }
  }
  ErrForR = rep(0, length(RSeq))
  RVec = rep(0, length(RSeq))
  for (r in 1:length(RSeq)) {
    RVec[r] = CVErrList[[r]]$R
    ErrForR[r] = min(CVErrList[[r]]$CVMat)
  }
  R = RVec[which.min(ErrForR)]
  CVScoreR = CVErrList[[which.min(ErrForR)]]$CVMat
  BestParaR = which(CVScoreR == min(CVScoreR), arr.ind = TRUE)
  eta1 = eta1Seq[BestParaR[1]]
  eta2 = eta2Seq[BestParaR[2]]
  return(list("CVScoreList" = CVErrList,"R" = R, "eta1" = eta1, "eta2" = eta2))
}


PCTraining = function(Theta0,R, eta1, eta2,SNeList,Omega){
  TrainPCPreList = TrainPCPrepare(Theta0 = Theta0,SNeList = SNeList)
  SNeList = TrainPCPreList$SNeList
  AtAList = TrainPCPreList$AtAList
  AtYList = TrainPCPreList$AtYList
  S = length(SNeList)
  PCRes = ADAMCpp(TrainIndex = 1:S, R = R, eta1 = eta1, eta2 = eta2,SNeList =SNeList, AtAList = AtAList,AtYList = AtYList,Omega = Omega)
  return(PCRes)
}

######################################


