# The script below estimates selection coefficients of of complex structural 
# variants (cxSV) and normal structural variants (SV)

##########################################
#                                        #
#     Load packages                      #
#                                        #
##########################################

# Load packages
library(GenomicRanges)
library(pracma)
library(RColorBrewer)

# Source functions
# AllFunctions <- list.files('/Users/heinrichzudohna/Library/Mobile Documents/com~apple~CloudDocs/cxSV/Functions/',
#                            full.names = T)
AllFunctions <- list.files('/home/hb54/cxSV/Functions/',
                           full.names = T)
sapply(AllFunctions, source)

##########################################
#                                        #
#     Set parameters                     #
#                                        #
##########################################

# False discovery rate for selected L1
FDR <- 0.1

# Human effective population size
PopSize <- 2*10^4

# Sample size to randomly sample standard SVs
SampleSizeSV <- 2000

# Specify file paths
# cxSVPath         <- '/Users/heinrichzudohna/Documents/cxSVData/8493cxSV_updatedinfo_AF_n_indv.txt'
# SVPath           <- '/Users/heinrichzudohna/Documents/cxSVData/simpleSV_combined_updatedinfo_AF_n_indv.txt'
# RegrOutputPath   <- paste0("/Users/heinrichzudohna/Library/Mobile Documents/com~apple~CloudDocs/cxSV/Results/SelectionRegressionResults_PopSize",
#                            PopSize, ".RData")
# AICTabOutputPath <- paste0("/Users/heinrichzudohna/Library/Mobile Documents/com~apple~CloudDocs/cxSV/Results/SelectionRegressionAICTab_PopSize",
#                            PopSize, ".csv")


cxSVPath         <- '/home/hb54/cxSVData/8493cxSV_updatedinfo_AF_n_indv.txt'
SVPath           <- '/home/hb54/cxSVData/simpleSV_combined_updatedinfo_AF_n_indv.txt'
RegrOutputPath   <- paste0("/home/hb54/cxSV/Results/SelectionRegressionResults_PopSize",
                           PopSize, ".RData")
AICTabOutputPath <- paste0("/home/hb54/cxSV/Results/SelectionRegressionAICTab_PopSize",
                           PopSize, ".csv")



##########################################
#                                        #
#     Load and process data              #
#                                        #
##########################################

cat("\n\nLoading and processing data ...")

# Read in files with cxSVs and standardSVs 
cxSVs <- read.delim(cxSVPath)
SVs   <- read.delim(SVPath)
SVSampleRows <- sample(1:nrow(SVs), SampleSizeSV)
SVSample <- SVs[SVSampleRows, ]
hist(cxSVs$AF, breaks = seq(0, 1, 0.001),
     xlim = c(0, 0.1), ylim = c(0, 1000))
hist(SVs$AF, breaks = seq(0, 1, 0.001),xlim = c(0, 0.1), ylim = c(0, 1000))

# Create a new file that combines both types of variants
Vars2Include <- c("span", "AF", "n_ind")
colnames(cxSVs)
AllSV <- rbind(cxSVs[,Vars2Include], SVSample[,Vars2Include])
AllSV$complex <- F
AllSV$complex[1:nrow(cxSVs)] <- T

# Get number of individuals ascertained
N <- AllSV$n_ind / AllSV$AF
N <- round(N, 0)
cat("done!\n\n")

##########################################################
#                                                        #
#   Fit effect of variant length and type on selection   #
#                                                        #
##########################################################

# Create a matrix of predictor variables (indicator for cxSV and size)
PredictMat <- AllSV[, c("complex", "span")]

cat("\n********   Estimating effect of variant length and type on selection   **********\n")
ModelFit_pracma <- FitSelectionModels_pracma(PredictMat,  
                                Freqs = AllSV$n_ind, 
                                Counts = rep(1, nrow(AllSV)), 
                                PopSize = PopSize, 
                                SampleSize = N,
                                aBorder = 0.003, 
                                bBorder = 10^(-1), 
                                cBorder = 10^(-6))


# Calculate the SV size difference, equivalent to the difference between SV and 
# cxSV
DeltaSize <- ModelFit_pracma$ML_abc$par["b.b"] / ModelFit_pracma$ML_abc$par["c.c"]

# Save everything
save.image(RegrOutputPath)
write.csv(ModelFit_pracma$AICTab, AICTabOutputPath)


##########################################################
#                                                        #
#   Fit different variant length effects on selection    #
#                                                        #
##########################################################

load(RegrOutputPath)

# Create a matrix of predictor variables (either size of cxSV or SV)
PredictMat <- cbind(AllSV$complex * AllSV$span, (1 - AllSV$complex) * AllSV$span)
colnames(PredictMat) <- c("SlopeComplex", "SlopeSimple")

cat("Fitting model with different slopes for different types of variants ...\n")
ML_DiffSlopes <- constrOptim(theta = c(a = ModelFit_pracma$ML_abc$par["a.a.a"], 
                                b = ModelFit_pracma$ML_abc$par["c.c"], 
                                c = ModelFit_pracma$ML_abc$par["c.c"]),
                      f = function(x) -AlleleFreqLogLik_3Par_pracma(
                        Freqs = AllSV$n_ind, 
                        Counts = rep(1, nrow(AllSV)), 
                        Predict = PredictMat,
                        a = x[1], b = x[2], c = x[3], N = PopSize, 
                        SampleSize = N),
                      grad = NULL,
                      ui = rbind(c(1, 0, 0),  c(0, 1, 0),  c(0, 0, 1), 
                                 c(-1, 0, 0), c(0, -1, 0), c(0, 0, -1)),
                      ci = c(a = -0.003, b = -10^(-6), c = -10^(-6), 
                             a = -0.003, b = -10^(-6), c = -10^(-6)),
                      method = "Nelder-Mead")
cat("done!\n\n")

##########################################################
#                                                        #
#      Fit logarithmic length effects on selection       #
#                                                        #
##########################################################

# Create a matrix of predictor variables (indicator for cxSV and logarithmic
# size)
PredictMat <- AllSV[, c("complex", "span")]
PredictMat$span <- log(PredictMat$span)

cat("Fitting model with different slopes for different types of variants ...\n")
ML_Log <- constrOptim(theta = c(a = ModelFit_pracma$ML_abc$par["a.a.a"], 
                                b = ModelFit_pracma$ML_abc$par["b.b"], 
                                c = ModelFit_pracma$ML_abc$par["c.c"]),
                             f = function(x) -AlleleFreqLogLik_3Par_pracma(
                               Freqs = AllSV$n_ind, 
                               Counts = rep(1, nrow(AllSV)), 
                               Predict = PredictMat,
                               a = x[1], b = x[2], c = x[3], N = PopSize, 
                               SampleSize = N),
                             grad = NULL,
                             ui = rbind(c(1, 0, 0),  c(0, 1, 0),  c(0, 0, 1), 
                                        c(-1, 0, 0), c(0, -1, 0), c(0, 0, -1)),
                             ci = c(a = -0.003, b = -10^(-2), c = -10^(-2), 
                                    a = -0.003, b = -10^(-2), c = -10^(-2)),
                             method = "Nelder-Mead")
cat("done!\n\n")

# Save everything
save.image(RegrOutputPath)

##########################################################
#                                                        #
#                    Plot results                        #
#                                                        #
##########################################################

# Best model parameters
Intercept_SV   <- ModelFit_pracma$ML_abc$par["a.a.a"] * PopSize
Intercept_cxSV <- sum(ModelFit_pracma$ML_abc$par[c("a.a.a", "b.b")]) * PopSize
Slope <- ModelFit_pracma$ML_abc$par["c.c"] * PopSize

# Compare size distribution of the two types of variants
BinWidth <- 100
MaxSize  <- 10000
Cols <- brewer.pal(3, "Dark2")
HistSV <- hist(SVs$span, breaks = seq(0, 1.1*max(SVs$span), BinWidth),
               xlim = c(0, 5*10^3))
HistcxSV <- hist(cxSVs$span, breaks = seq(0, 1.1*max(SVs$span), BinWidth),
                 xlim = c(0, 5*10^3))
pdf(file = "/home/hb54/cxSV/Figures/Selection coefficient vs variant size.pdf")
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(c(0, MaxSize), c(0, max(c(HistSV$density, HistcxSV$density))), 
     xlab = "Variant size [bp]", ylab  = "Density",
     type = "n")
rect(xleft = c(0, HistSV$mids[-length(HistSV$mids)] + 0.5 * BinWidth), 
     ybottom = 0, xright = c(HistSV$mids) - 0.1, ytop = HistSV$density,
     col = Cols[1])
rect(xleft = HistcxSV$mids[-length(HistcxSV$mids)], 
     ybottom = 0, xright = HistcxSV$mids + 0.5 * BinWidth, 
     ytop = HistcxSV$density,
     col = Cols[2])

# Add lines of fitness values
par(new = TRUE)
xVals <- 0:MaxSize
MinFit <- Intercept_cxSV + Slope * MaxSize
plot(xVals, Intercept_SV + Slope * xVals, type = "l",              
     axes = FALSE, xlab = "", ylab = "", ylim = c(MinFit, 0), col = Cols[1])
lines(xVals, Intercept_cxSV + Slope * xVals, col = Cols[2])
axis(side = 4, at = pretty(c(MinFit, 1)))      
mtext(expression(N[e]~"*"~s), side = 4, line = 3) 

# Add legend
legend("topright", bty = "n", col = Cols[1:2], legend = c("SV", "cxSV"), lty = 1)
dev.off()