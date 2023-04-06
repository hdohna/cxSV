# The following script generates and runs slurm scripts to estimate selection
# coefficient models for SVs and cxSV

# Source functions
AllFunctions <- list.files('/home/hb54/cxSV/Functions/', full.names = T)
sapply(AllFunctions, source)

# Command to load r
cmdLoadR  <- 'module load R'
cmdScript <- 'Rscript /home/hb54/cxSV/Scripts/EstimateSelectionPars_InputFile.R'

# Function to set up and submita run
SubmitRun <- function(PopSize, SVPath, Cols2Include, ThetaStart){
  
  # Generate time stamp for input file
  InputPath   <- paste0("/home/hb54/cxSV/Data/In_PopSize",
                             PopSize,"_Cols_", paste(Cols2Include, collapse = "_"),
                             ".RData")
  
  # Save parameters into input file
  save(list = c("PopSize", "SVPath", "Cols2Include", "ThetaStart"), 
       file = InputPath)
  
  # Create and run script
  ScriptFile <- paste0('/home/hb54/cxSV/SelScript_Pop', PopSize,"_Cols_", 
                       paste(Cols2Include, collapse = "_"))
  Cmds <- c(cmdLoadR, paste(cmdScript, InputPath))
  CreateAndCallSlurmScript(file = ScriptFile,
                           RunTime = '24:00:00',
                           Mem = '20G',
                           SlurmCommandLines =  Cmds)
  
}

# Submit runs for individual models that did not complete

# 1-parameter model for population size of 20,000
SubmitRun(PopSize = 2*10^4, 
          SVPath = "/home/hb54/cxSVData/sampledSV_Tue_Mar_28_08:33:50_2023.csv",
          Cols2Include = 0,
          ThetaStart = c(a = -0.0004))

# linear 3-parameter model for population size of 20,000
SubmitRun(PopSize = 2*10^4, 
          SVPath = "/home/hb54/cxSVData/sampledSV_Tue_Mar_28_08:33:50_2023.csv",
          Cols2Include = 1:2,
          ThetaStart = c(a = -0.0003921515, b = -0.0002227735, c =  -9.067496e-08))

# linear 2-slope model for population size of 20,000
SubmitRun(PopSize = 2*10^4, 
          SVPath = "/home/hb54/cxSVData/sampledSV_Tue_Mar_28_08:33:50_2023.csv",
          Cols2Include = 4:5,
          ThetaStart = c(a = -0.002144079, b = -3.404487e-07, c =  -3.0895e-07))

# log-linear 3-parameter model for population size of 20,000
SubmitRun(PopSize = 2*10^4, 
          SVPath = "/home/hb54/cxSVData/sampledSV_Tue_Mar_28_08:33:50_2023.csv",
          Cols2Include = c(1, 3),
          ThetaStart = c(a = -4.979167e-05, b = -0.0001997917, c =  5.25e-07))



