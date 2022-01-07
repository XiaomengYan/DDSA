# Output Functons
Mean_SEDTraining = function(SNeList, Omega,
                            lambda1 = 1e-5, lambda2 = 1e-3,
                            isCV = FALSE, KFold = 5, lambda1Seq = exp(seq(-10,3,length.out = 10)),
                            lambda2Seq = exp(seq(-10,3,length.out = 10)), outPath = "CVScoreMat_Mean.dat",
                            isPara = TRUE, NoCore = 5){
  if(isCV){
    CVMeanResList = CV.MeanTraining(lambda1Seq = lambda1Seq,lambda2Seq = lambda2Seq,SNeList = SNeList,Omega = Omega,
                                    KFold = KFold, ifParallel = isPara, NoCore = NoCore)
    CVScoreMat = CVMeanResList$CVScoreMat
    write.table(CVScoreMat, file = outPath ,quote = FALSE, row.names = FALSE, col.names = FALSE)
    lambda1 = CVMeanResList$lambda1
    lambda2 = CVMeanResList$lambda2
  }
  cat("Begin to Train Mean Surface with Best Parameters!")
  Theta0 = MeanTraining(lambda1 = lambda1,lambda2 = lambda2,SNeList = SNeList,Omega = Omega)
  return(list("Theta0"  = Theta0, "paras" = list("lambda1" = lambda1, "lambda2" = lambda2)))
}

PC_SEDTraining = function(SNeList, Omega, Theta0,
                          R = 2, eta1 = 1e-6, eta2 = 1e-3,
                          isCV = FALSE,KFold = 5, RSeq = 1:20, eta1Seq = exp(seq(-10,3,length.out = 10)),
                          eta2Seq = exp(seq(-10,3,length.out = 10)), outPath = "CVPCScore.RData",
                          isPara = TRUE, NoCore = 5){
  if(isCV){
    CVPCResList = CV.PCTraining(Theta0 = Theta0,RSeq = RSeq, eta1Seq = eta1Seq,eta2Seq = eta2Seq,SNeList = SNeList,Omega = Omega,
                                KFold = KFold,ifParallel = isPara,NoCore = NoCore)
    CVScoreList = CVPCResList$CVScoreList
    save(CVScoreList ,file = outPath)
    R = CVPCResList$R
    eta1 = CVPCResList$eta1
    eta2 = CVPCResList$eta2
  }
  PCRes = PCTraining(Theta0 = Theta0,R = R, eta1 = eta1,eta2 = eta2,SNeList = SNeList,Omega = Omega)
  return(list("pcResult" = PCRes, "paras" = list("R" = R, "eta1" = eta1, "eta2" = eta2)))
}
