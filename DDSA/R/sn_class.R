# supernova class
sn_class <- R6Class(
  "sn_class",
  portable = FALSE,
  public = list(
    Name = NA,
    N = numeric(0),
    PhaseVec = numeric(0),
    FluxList = numeric(0),
    FluxMean0List = numeric(0),
    WavelengthList = numeric(0),
    ErrList = numeric(0),
    WmatList = numeric(0),
    BmatList = numeric(0),

    initialize = function(SNedf, SNe,WaveletObj, SplineObj){
      Name <<- SNe
      PhaseVec <<- unique(SNedf$Phase)
      N <<- length(PhaseVec)
      SplitList  = split(SNedf, f = SNedf$Phase)
      FluxList = list()
      WavelengthList = list()
      ErrList = list()
      WmatList = list()
      BmatList = list()

      for (n in 1:N) {
        FluxList[[n]] = SplitList[[n]]$Flux
        WavelengthList[[n]] = SplitList[[n]]$Wavelength
        ErrList[[n]] = SplitList[[n]]$Err
        WmatList[[n]] = WaveletObj$WaveletEvalAt(WavelengthList[[n]])
        BmatList[[n]] = matrix(SplineObj$bD0EvalAt(PhaseVec[n]),ncol = 1)
      }

      FluxList <<- FluxList
      WavelengthList <<- WavelengthList
      ErrList <<- ErrList
      WmatList <<- WmatList
      BmatList <<- BmatList
    },

    RemoveMean = function(Theta0){
      FluxMean0List = list()
      for (n in 1:N) {
        FluxMean0List[[n]] = FluxList[[n]] - WmatList[[n]]%*%Theta0%*%BmatList[[n]]
      }
      FluxMean0List <<- FluxMean0List
    }

  ),
  private = list()
)


# LogSNeClass <- R6Class(
#   "LogSNeClass",
#   portable = FALSE,
#   public = list(
#     Name = NA,
#     N = numeric(0),
#     PhaseVec = numeric(0),
#     FluxList = numeric(0),
#     FluxMean0List = numeric(0),
#     WavelengthList = numeric(0),
#     ErrList = numeric(0),
#     WmatList = numeric(0),
#     BmatList = numeric(0),
#
#     initialize = function(SNedf, SNe,WaveletObj, SplineObj){
#       Name <<- SNe
#       PhaseVec <<- unique(SNedf$Phase)
#       N <<- length(PhaseVec)
#       SplitList  = split(SNedf, f = SNedf$Phase)
#       FluxList = list()
#       WavelengthList = list()
#       ErrList = list()
#       WmatList = list()
#       BmatList = list()
#
#       for (n in 1:N) {
#         FluxList[[n]] = log(SplitList[[n]]$Flux+1)
#         WavelengthList[[n]] = SplitList[[n]]$Wavelength
#         ErrList[[n]] = SplitList[[n]]$Err/abs(SplitList[[n]]$Flux+1)
#         WmatList[[n]] = WaveletObj$WaveletEvalAt(WavelengthList[[n]])
#         BmatList[[n]] = matrix(SplineObj$bD0EvalAt(PhaseVec[n]),ncol = 1)
#       }
#
#       FluxList <<- FluxList
#       WavelengthList <<- WavelengthList
#       ErrList <<- ErrList
#       WmatList <<- WmatList
#       BmatList <<- BmatList
#     },
#
#     RemoveMean = function(Theta0){
#       FluxMean0List = list()
#       for (n in 1:N) {
#         FluxMean0List[[n]] = FluxList[[n]] - WmatList[[n]]%*%Theta0%*%BmatList[[n]]
#       }
#       FluxMean0List <<- FluxMean0List
#     }
#
#   ),
#   private = list()
# )
