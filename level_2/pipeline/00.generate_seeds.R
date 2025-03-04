args <- commandArgs(TRUE)

nseeds <- as.numeric(args[1])
globalseed <- as.numeric(args[2])

set.seed(globalseed)

writeLines(text=as.character(sample.int(1e6, nseeds)),
           con="seedfile")
