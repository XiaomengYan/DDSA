
# Orthonormal basis function, its first and second derivative
# basisD0_fun, basisD1_fun, basisD2_fun
# Omega matrix for smoothing
R_spline_class <- R6Class(
  "R_spline_class",
  portable = FALSE,
  public = list(
    phase_min = NA,
    phase_max = NA,
    Omega = numeric(0),
    bknots = NA,

    initialize = function(pmin, pmax){
      phase_min <<- pmin
      phase_max <<- pmax
      createBFList()
    },

    bD0EvalAt = function(phase){
      fChoice = basisD0FList
      basisRes = bDXEvalAt(phase, fChoice)
      return(basisRes)
    },

    bD1EvalAt = function(phase){
      fChoice = basisD1FList
      basisRes = bDXEvalAt(phase, fChoice)
      return(basisRes)
    },

    bD2EvalAt = function(phase){
      fChoice = basisD2FList
      basisRes = bDXEvalAt(phase, fChoice)
      return(basisRes)
    }

  ),## public

  private = list(
    nTXX = 1e3,
    basisD0FList = list(),
    basisD1FList = list(),
    basisD2FList = list(),

    bDXEvalAt = function(phase, fChoice){
      phase4basis = (phase - phase_min) / (phase_max - phase_min)
      numBasis = length(fChoice)
      basisRes = sapply(1:numBasis,
                        FUN = function(j) fChoice[[j]](phase4basis))
      return(basisRes)
    },


    createBFList = function(){
      bSL = genBSplineList()
      oSL = genOSplineList(bSL)
      computeOmega(oSL)
      genSFunList(oSL)
    },

    computeOmega = function(oSL){
      stepSize = (phase_max - phase_min)/nTXX
      OmegaCC = t(oSL[[3]]) %*% oSL[[3]]* stepSize
      Omega <<- OmegaCC/norm(OmegaCC, "F")
    },

    genBSplineList = function(){
      # Version 1
      # bknots <<- c(c(-0.1,-0.05),
      #              seq(0, 1, length.out = 15),
      #              c(1.05,1.1))
      # Version 2
      # bknots <<- c(c(-0.075,-0.05,-0.025),
      #              seq(0, 1, length.out = 21),
      #              c(1.025,1.05,1.075))

      # Version 3
      bknots <<- c(c(-0.375,-0.25,-0.125),
                   seq(0, 1, length.out = 5),
                   c(1.125,1.25,1.375))

      tXX = seq(0, 1, length.out = nTXX)
      bSplineXList = list()
      for (j in 1:3){
        bSplineXList[[j]] =
          splineDesign(knots = bknots,
                       x = tXX, ord = 4,
                       derivs = rep(j - 1, nTXX))
      }
      return(bSplineXList)
    },


    genOSplineList = function(bSplineXList){
      ## Orthogonalize B-spline
      oSplineXList = list()
      R = qr.R( qr(bSplineXList[[1]]) )
      RInv = solve(R)
      stepSize = (phase_max - phase_min) * 1/nTXX
      oConst = sqrt(1/stepSize)
      ## /int b^2(t) dt = \sum b^2(t_k) * stepSize = 1 (normalization)
      ## b(t_k) =
      for (j in 1:3){
        oSplineXList[[j]] = bSplineXList[[j]] %*% RInv* oConst
      }
      return(oSplineXList)
    },

    genSFunList = function(oSL){
      ## Convert Basis (d0 d1 d2) to functions
      dimBasis = dim(oSL[[1]])[2]
      tXX = seq(0, 1, length.out = nTXX)
      for (j in 1:dimBasis){
        basisD0FList[[j]] <<- approxfun(x = tXX,
                                        y = oSL[[1]][,j])
        basisD1FList[[j]] <<- approxfun(x = tXX,
                                        y = oSL[[2]][,j])
        basisD2FList[[j]] <<- approxfun(x = tXX,
                                        y = oSL[[3]][,j])
      }
    }

  )## private
)## basis Class


