Args <- commandArgs(trailingOnly = T)
cat("Input argument is", eval(parse(text = Args[1])))
