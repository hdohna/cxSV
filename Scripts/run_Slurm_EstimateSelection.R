# The following script generates and runs slurm scripts to estimate selection
# coefficient models for SVs and cxSV

# Source functions
AllFunctions <- list.files('/home/hb54/cxSV/Functions/', full.names = T)
sapply(AllFunctions, source)

# Command to load r
cmdLoadR <- 'module load R'
cmdScript <- 'Rscript /home/hb54/cxSV/Scripts/EstimateSelectionPars_Module.R'

# Specify number of SVs to sample
NSampleSV <- 2000
PopSize <- 5*10^3
for (PopSize in c(5, 20) * 10^3){
  
  # Model 1: a single selection coefficient
  ScriptFile <- '/home/hb54/cxSV/Functions/SelMod1'
  Cmds <- c(cmdLoadR, paste(cmdScript, PopSize, NSampleSV, 0))
  CreateAndCallSlurmScript(file = ScriptFile,
                           RunTime = '24:00:00',
                           Mem = '20G', 
                           SlurmCommandLines =  Cmds)
                           
  # Model 2: selection coefficient differs between cxSV and SV
  ScriptFile <- '/home/hb54/cxSV/Functions/SelMod2'
  Cmds <- c(cmdLoadR, paste(cmdScript, PopSize, NSampleSV, 1))
  CreateAndCallSlurmScript(file = ScriptFile,
                           RunTime = '24:00:00',
                           Mem = '20G', 
                           SlurmCommandLines =  Cmds)
  
  # Model 3: selection coefficient depends on size
  ScriptFile <- '/home/hb54/cxSV/Functions/SelMod3'
  Cmds <- c(cmdLoadR, paste(cmdScript, PopSize, NSampleSV, 2))
  CreateAndCallSlurmScript(file = ScriptFile,
                           RunTime = '24:00:00',
                           Mem = '20G', 
                           SlurmCommandLines =  Cmds)
  
  # Model 4: selection coefficient differs between cxSV and SV and depends on size
  ScriptFile <- '/home/hb54/cxSV/Functions/SelMod4'
  Cmds <- c(cmdLoadR, paste(cmdScript, PopSize, NSampleSV, '"c(1,2)"'))
  CreateAndCallSlurmScript(file = ScriptFile,
                           RunTime = '24:00:00',
                           Mem = '20G', 
                           SlurmCommandLines =  Cmds)

  # Model 5: selection coefficient differs between cxSV and SV and depends on 
  # logarithmic size
  ScriptFile <- '/home/hb54/cxSV/Functions/SelMod5'
  Cmds <- c(cmdLoadR, paste(cmdScript, PopSize, NSampleSV, '"c(1,3)"'))
  CreateAndCallSlurmScript(file = ScriptFile,
                           RunTime = '24:00:00',
                           Mem = '20G', 
                           SlurmCommandLines =  Cmds)

  # Model 6: size dependence of selection coefficient differs between cxSV and SV
  ScriptFile <- '/home/hb54/cxSV/Functions/SelMod6'
  Cmds <- c(cmdLoadR, paste(cmdScript, PopSize, NSampleSV, '"c(4,5)"'))
  CreateAndCallSlurmScript(file = ScriptFile,
                           RunTime = '24:00:00',
                           Mem = '20G', 
                           SlurmCommandLines =  Cmds)
  
}


