#### Wavelet Basis Function

R_wavelet_class <- R6Class(
  "R_wavelet_class",
  portable = FALSE,
  public = list(
    lambda_min = NA,
    lambda_max = NA,
    NumBasis = NA,

    # Initialization (Input: lambda_min, lambda_max)
    initialize = function(lambdamin,lambdamax){
      lambda_min <<- lambdamin
      lambda_max <<- lambdamax
      creatWFlist()
    },
    WaveletEvalAt = function(lambdaSeq){
      basisRes = WFEvalAt(lambdaSeq)
      return(basisRes)
    }
  ),## public

  private = list(
    WFlist = list(),

    WFEvalAt = function(lambdaSeq){
      for (p in 1:NumBasis) {
        res = matrix(0,ncol = NumBasis,nrow = length(lambdaSeq))
        for (p in 1:NumBasis) {
          res[,p] = WFlist[[p]](lambdaSeq)
        }
        return(res)
      }
    },

    creatWFlist = function(){
      #library(SNeIaSurf)
      waveletObj = new(wavelet_class)
      #lambda_min = 3500
      #lambda_max = 9000
      WFlist = list()
      NumBasis <<-  waveletObj$set_parameter(lambda_min,lambda_max,500,4)
      LambdaSeq = seq(from = lambda_min,to = lambda_max, by = 1)
      Wmat = waveletObj$compute_wavelet_mat(LambdaSeq)
      for (p in 1:NumBasis) {
        WFlist[[p]] = approxfun(x = LambdaSeq,y = Wmat[,p])
      }
      WFlist <<- WFlist
    }

  )# private
)# Wavelet Class
