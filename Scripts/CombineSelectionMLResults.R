# The following script reads in results of individual runs that estimated 
# maximum likelihood parameters of selection coefficient models

# Get all files with results:
ResultFiles <- list.files("cxSV/Results/", 
                          pattern = "SelectionRegressionResults_", full.names = T)

# Function to find first non-zero digit behind the decimal point in numbers < 1
FirstNonZero <- function(x, MaxDigits  = 10){
  min(which((abs(x) - 10^-(0:MaxDigits)) > 0))
}

# Function to process results
ResultPath <- ResultFiles[1]
ResultPath
load(ResultPath)
idxVars
ML$par
Args[3]
Cols2Include
ResultTable <- t(sapply(ResultFiles[-1], function(ResultPath){
  load(ResultPath)
  # Correct a mistake in the original script EstimateSelectionPars_Module.R here
  idxVars <- ifelse(length(idxVars) == 1, 1, c(1, idxVars + 1))
  VarNames <- c("Int", "Complex", "Size", "LogSize", "SlopeCx", "SlopeSV")[idxVars]
  ParsSimple <- sapply(ML$par, function(x){
    NDigit <- FirstNonZero(x)
    round(x, digits = NDigit)
  })
  c(
    Likelihood = round(ML$value, 2),
    FittedParameters = paste(VarNames, collapse = ", "),
    NrParameters = length(ML$par),
    AIC = round(2 * (length(ML$par) + ML$value), 2),
    Results = paste(VarNames, ParsSimple, sep = " = ",
                    collapse = ", ")
    
  )
}))
ResultTable

