##############################################
#
# General description:
#
#   The following function calculates the log probability of the count value
#   for an allele under selection (equation after Boissinot at al. 2006 PNAS)

# Input:
#
#     k: number of samples carrying the allele
#     s: selection coefficient
#     N: population size
#     SampleSize: sample size

# Comment:
#     This function requires the function AlleleFreqTime and package pracma

##############################################

AlleleFreqSample_pracma <- function(k, s, N, SampleSize = 2504, DetectProb = 1){
    

    # Calculate integration constant
    IntConst <- integral(function(x) {
        AlleleFreqTime(x, s, N) 
      },  0, 1) - 
    integral(function(x) {
        AlleleFreqTime(x, s, N) * (1 - DetectProb * x)^SampleSize 
      },  0, 1)
    
    # Calculate probability of obtaining k alleles in a sample of size 
    # SampleSize
    log(
      integral(function(x) {
          AlleleFreqTime(x, s, N) * 
          dbinom(k, SampleSize, DetectProb * x)
        }, 0, 1) - 
      integral(function(x) {
          AlleleFreqTime(x, s, N) * (1 - DetectProb * x)^SampleSize *  
            dbinom(k, SampleSize, DetectProb * x)
        }, 0, 1)
      ) -
      log(IntConst)
    # The lines below are for deletions relative to the reference

  
}


