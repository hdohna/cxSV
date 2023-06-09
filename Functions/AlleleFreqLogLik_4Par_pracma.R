##############################################
#
# General description:
#
#   The following function calculates the probability log-likelihood of sampled
#   frequencies of alleles under selection (equation after Boissinot at al. 2006 PNAS)

# Input:
#
#     Freqs:  vector of observed allele frequencies (integer)
#     Counts: vector of counts (number of alleles with the frequencies 
#             specified in Freqs)
#     Predict: matrix (n x 2) of predictor variable values
#     a:      intercept
#     b:      slope for first predictor
#     c:      slope for second predictor
#     N:      population size

# Comment:
#     This function requires the function AlleleFreqTime and AlleleFreqSample

##############################################

AlleleFreqLogLik_4Par_pracma <- function(Freqs, Counts, Predict, a, b, c, d, SD = NULL, N, 
                                  SampleSize = 2504,
                                  MinFactor = 2, 
                                  NrObsExtrapol = 10, 
                                  verbose = T, showInfIndices = F,
                                  LowerS = -1,
                                  UpperS = 1){
  if ((length(Freqs) != length(Counts)) | (length(Freqs) != nrow(Predict)) |
      (nrow(Predict) != length(Counts)) ){
    browser()
    stop("Freqs, Counts and Predict vector have to have the same length\n")
  }
  if (length(SampleSize) == 1){
    SampleSize <- rep(SampleSize, length(Freqs))
  }
  sVals <- as.matrix(Predict) %*% c(b, c, d)
  LogLikVals <- sapply(1:length(Freqs), function(i){
      Counts[i] * AlleleFreqSample_pracma(Freqs[i], a + sVals[i], N, 
                                   SampleSize = SampleSize[i])
      
  })

  blnInf       <- is.infinite(LogLikVals)
#  LogLikVals[blnInf] <- MinFactor*min(LogLikVals[!blnInf])

  if (verbose){
    cat("parameter values: a =", a, "b =", b, "c = ", c, "d = ", d,
        "SD = ", SD,
        "log-likelihood =", sum(LogLikVals), "\n")
    cat("Number of observations with infinite log-likelihood:", sum(blnInf), "\n")
  }
  if (showInfIndices){
    cat("Indices of infinte values:", paste(which(blnInf), collapse = "\n"), "\n")
  }
  sum(LogLikVals)
}


