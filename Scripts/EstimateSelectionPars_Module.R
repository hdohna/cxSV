# The script below estimates selection coefficients of of complex structural 
# variants (cxSV) and normal structural variants (SV)

# Collect arguments:
#    first argument: effective population size
#    second argument: sample size of standard SVs
#    third argument: columns of predictor matrix to include (0 means no column)
Args <- commandArgs(trailingOnly = T)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

cat("Setting parameters ... ")
# Human effective population size
PopSize <- as.numeric(Args[1])

# Sample size to randomly sample standard SVs
SampleSizeSV <- as.numeric(Args[2])

# Get the columns of the prediction matrix to include
Cols2Include <- eval(parse(text = Args))
if(!all(Cols2Include %in% 0:2)) stop("Third argument should be integers between 0 and 5\n")

# Define borders for constrained parameter estimations
aBorder = 0.003 
bBorder = 10^(-1) 
cBorder = 10^(-6)

# Specify file paths
FunctionPath     <- '/home/hb54/cxSV/Functions/'
cxSVPath         <- '/home/hb54/cxSVData/8493cxSV_updatedinfo_AF_n_indv.txt'
SVPath           <- '/home/hb54/cxSVData/simpleSV_combined_updatedinfo_AF_n_indv.txt'
idxVars <- ifelse(Cols2Include == 0, 1, c(1, Cols2Include))
VarNames <- c("Int", "Complex", "Size", "LogSize")[idxVars]
RegrOutputPath   <- paste0("/home/hb54/cxSV/Results/SelectionRegressionResults_PopSize",
                           PopSize,"_Vars_", paste(VarNames, collapse = "_"),
                          ".RData")
cat(" done!\n") 
cat("Model will fit coefficient for:", VarNames, "\n")

##########################################
#                                        #
#     Load packages and functions        #
#                                        #
##########################################

cat("Loading packages and functions ... ")
# Load packages
library(pracma)

# Source functions
AllFunctions <- list.files(FunctionPath, full.names = T)
sapply(AllFunctions, source)
cat(" done!\n") 

##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("\n\nLoading and processing data ...")

# Read in files with cxSVs and standardSVs 
cxSVs        <- read.delim(cxSVPath)
SVs          <- read.delim(SVPath)
SVSampleRows <- sample(1:nrow(SVs), SampleSizeSV)
SVSample     <- SVs[SVSampleRows, ]

# Create a new file that combines both types of variants
Vars2Include <- c("span", "AF", "n_ind")
AllSV        <- rbind(cxSVs[,Vars2Include], 
                      SVSample[,Vars2Include])

# Add an indicator variable whether SV is complex 
AllSV$complex <- F
AllSV$complex[1:nrow(cxSVs)] <- T

# Add logarithmic span
AllSV$logspan <- log(AllSV$span)

# Get number of individuals ascertained
NSample <- AllSV$n_ind / AllSV$AF
NSample <- round(NSample, 0)
cat("done!\n\n")

##########################################################
#                                                        #
#   Fit effect of variant length and type on selection   #
#                                                        #
##########################################################

# Create a matrix of predictor variables (indicator for cxSV and size)
PredictMat <- cbind(AllSV[, c("complex", "span", "logspan")],
                    AllSV$complex * AllSV$span, 
                    (1 - AllSV$complex) * AllSV$span)
colnames(PredictMat) <- c("complex", "span", "logspan", "slope_complex",
                          "slope_simple")
if(Cols2Include == 0){
  cat("****. Estimating maximum likelihood for a single selection coefficient  ****\n")
  ML <-  constrOptim(theta = c(a = 0),
                       f = function(x) -AlleleFreqLogLik_3Par_pracma(
                         Freqs = AllSV$n_ind,
                         Counts = rep(1, nrow(AllSV)),
                         Predict = PredictMat,
                         a = x[1], 
                         b = 0, 
                         c = 0, 
                         N = PopSize,
                         SampleSize = NSample),
                       grad = NULL,
                       ui = rbind(1, -1),
                       ci = c(a = -aBorder, a = -aBorder),
                       method = "Nelder-Mead")
  cat("done!\n")
  
} else {
  PredictMat <- PredictMat[,Cols2Include]
  cat("Estimate effect of", colnames(PredictMat), "on selection ...\n")
  NPar <- ncol(PredictMat) + 1
  ThetaStart <- c(a = 0, b = 0, c = 0)[1:NPar]
  UI <- rbind(diag(3), -diag(3)) [c(1:NPar, 3 + 1:NPar),1:NPar]
  Border1 <- c(a = -aBorder, b = -bBorder, c = -cBorder)[c(1, 1+Cols2Include)]
  CI = c(Border1, Border1)
  
  ML <-  constrOptim(theta = c(a = ML_a$par, b = 0),
                        f = function(x) -AlleleFreqLogLik_3Par_pracma(
                          Freqs = Freqs,
                          Counts = rep(1, length(Freqs)),
                          Predict = PredictMat,
                          a = x[1], 
                          b = x[2], 
                          c = ifelse(NPar == 3, x[3], 0), 
                          N = PopSize, 
                          SampleSize = SampleSize),
                        grad = NULL,
                        ui = UI,
                        ci = CI,
                        method = "Nelder-Mead")
  cat("done!\n")
  
}

# Save everything
cat("Saving results ... ")
save.image(RegrOutputPath)
cat("done!")
