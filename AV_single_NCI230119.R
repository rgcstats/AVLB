source("AV_preamble_220119.R")
k <- commandArgs() ; k <- k[length(k)] ; k <- as.numeric(k)
sink(file=paste0("test",k,".rout"))
t1 <- proc.time()[3]
simres <- simulate.fn(k=k,Rdiv=1,scenarios=scenarios)
save(simres,file=paste0("simres",k,".rdata"))
t2 <- proc.time()[3]
t2-t1
