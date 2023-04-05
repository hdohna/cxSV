# The following script generates and runs slurm scripts to estimate selection
# coefficient models for SVs and cxSV

# Source functions
AllFunctions <- list.files('/home/hb54/cxSV/Functions/', full.names = T)
sapply(AllFunctions, source)

# Command to load r
cmdLoadR  <- 'module load R'
cmdScript <- 'Rscript /home/hb54/cxSV/Scripts/EstimateSelectionPars_Module.R'

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


