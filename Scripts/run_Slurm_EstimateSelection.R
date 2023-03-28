# The following script generates and runs slurm scripts to estimate selection
# coefficient models for SVs and cxSV

# Source functions
AllFunctions <- list.files('/home/hb54/cxSV/Functions/', full.names = T)
sapply(AllFunctions, source)

# Specify number of SVs to sample
NSampleSV <- 2000

# Read in SV data, draw a subset and random and write it out
SVs          <- read.delim('/home/hb54/cxSVData/simpleSV_combined_updatedinfo_AF_n_indv.txt')
SVSampleRows <- sample(1:nrow(SVs), NSampleSV)
SVSample     <- SVs[SVSampleRows, ]
TStamp <- timestamp(prefix = "", suffix = "")
TStamp <- gsub(" ", "_", TStamp)
SVOutPath    <- paste0('/home/hb54/cxSVData/sampledSV_',TStamp,".csv")
write.csv(SVSample, SVOutPath, row.names = F)

# Command to load r
cmdLoadR  <- 'module load R'
cmdScript <- 'Rscript /home/hb54/cxSV/Scripts/EstimateSelectionPars_Module.R'

PopSize <- 5*10^3
for (PopSize in c(5, 20) * 10^3){
  
  # Model 1: a single selection coefficient
  ScriptFile <- paste0('/home/hb54/cxSV/SelMod1_Pop', PopSize)
  Cmds <- c(cmdLoadR, paste(cmdScript, PopSize, SVOutPath, 0))
  CreateAndCallSlurmScript(file = ScriptFile,
                           RunTime = '24:00:00',
                           Mem = '20G',
                           SlurmCommandLines =  Cmds)

  # Model 2: selection coefficient differs between cxSV and SV
  ScriptFile <- paste0('/home/hb54/cxSV/SelMod2_Pop', PopSize)
  Cmds <- c(cmdLoadR, paste(cmdScript, PopSize, SVOutPath, 1))
  CreateAndCallSlurmScript(file = ScriptFile,
                           RunTime = '24:00:00',
                           Mem = '20G',
                           SlurmCommandLines =  Cmds)

  # Model 3: selection coefficient depends on size
  ScriptFile <- paste0('/home/hb54/cxSV/SelMod3_Pop', PopSize)
  Cmds <- c(cmdLoadR, paste(cmdScript, PopSize, SVOutPath, 2))
  CreateAndCallSlurmScript(file = ScriptFile,
                           RunTime = '24:00:00',
                           Mem = '20G',
                           SlurmCommandLines =  Cmds)

  # Model 4: selection coefficient differs between cxSV and SV and depends on size
  ScriptFile <- paste0('/home/hb54/cxSV/SelMod4_Pop', PopSize)
  Cmds <- c(cmdLoadR, paste(cmdScript, PopSize, SVOutPath, '"c(1,2)"'))
  CreateAndCallSlurmScript(file = ScriptFile,
                           RunTime = '24:00:00',
                           Mem = '20G',
                           SlurmCommandLines =  Cmds)

  # Model 5: selection coefficient differs between cxSV and SV and depends on 
  # logarithmic size
  ScriptFile <- paste0('/home/hb54/cxSV/SelMod5_Pop', PopSize)
  Cmds <- c(cmdLoadR, paste(cmdScript, PopSize, SVOutPath, '"c(1,3)"'))
  CreateAndCallSlurmScript(file = ScriptFile,
                           RunTime = '24:00:00',
                           Mem = '20G', 
                           SlurmCommandLines =  Cmds)

  # Model 6: size dependence of selection coefficient differs between cxSV and SV
  ScriptFile <- paste0('/home/hb54/cxSV/SelMod6_Pop', PopSize)
  Cmds <- c(cmdLoadR, paste(cmdScript, PopSize, SVOutPath, '"c(4,5)"'))
  CreateAndCallSlurmScript(file = ScriptFile,
                           RunTime = '24:00:00',
                           Mem = '20G', 
                           SlurmCommandLines =  Cmds)
  
}


