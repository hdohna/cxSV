# The following script reads in results of individual runs that estimated 
# maximum likelihood parameters of selection coefficient models

# Get all files with results:
ResultFiles <- list.files("cxSV/Results/", 
                          pattern = "SelectionRegressionResults_", full.names = T)

# Function to find first non-zero digit behind the decimal point in numbers < 1
FirstNonZero <- function(x, MaxDigits  = 10){
  min(which((abs(x) - 10^-(0:MaxDigits)) > 0))
}

# Process results from different runs and put them in a table
ResultTable <- t(sapply(ResultFiles[-1], function(ResultPath){
  load(ResultPath)
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

# Turn table into a data frame and make numeric columns numeric
ResultTable <- as.data.frame(ResultTable)
for (Col in c("Likelihood", "NrParameters", "AIC")){
  ResultTable[, Col] <- as.numeric(ResultTable[, Col])
}

write(ResultTable, "cxSV/Results/ResultTable.csv")

# Likelihood ratio test to compare top two models
orderAIC <- order(ResultTable$AIC)
DiffLik <- ResultTable$Likelihood[orderAIC[2]] - ResultTable$Likelihood[orderAIC[1]]
DiffPar   <- ResultTable$NrParameters[orderAIC[1]] - ResultTable$NrParameters[orderAIC[2]]
pchisq(q = 2*DiffLik, df = DiffPar, ncp = 0, lower.tail = F, log.p = FALSE)

