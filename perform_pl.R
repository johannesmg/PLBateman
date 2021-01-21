library(GenSA)
library(parallel)
library(gsl)
library(lhs)                                       
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(numDeriv)
library(minpack.lm)


source("functions.R")
source("gradient_functions.R")

args = commandArgs(trailingOnly=TRUE)

data.file=args[1]
in.file=args[2]
out.file=args[3]
parameter.bounds.file=args[4]
pl.par.name=args[5]
step.size.multiplier=as.numeric(as.character(args[6]))
no.cores=as.numeric(as.character(args[7]))
optim.method=args[8]
fnscale=as.numeric(as.character(args[9]))
factr=as.numeric(as.character(args[10]))
maxit=as.numeric(as.character(args[11]))


par.names <- c("g1p","g2p","k01","n1","k02","n2","a","b","g1","g2","n","dt")
no.pars=length(par.names)
parm.bounds <- load.bounds.from.file(parameter.bounds.file)
names(parm.bounds) <- c("lo","hi","start")
parm.bounds$start <- NA

min.cost=read_csv(in.file)
current.gene=min.cost$gene[1]

parm.bounds$start=as.numeric(min.cost[,par.names])
names(parm.bounds$start)=par.names

all.data=provide.data(data.file,current.gene)

if (optim.method=="lbfgsb"){
    optim.arg=list(trace=0,fnscale=fnscale,factr=factr,maxit=maxit)
}else if (optim.method=="lm"){
    optim.arg=nls.lm.control(ftol=1e-25,ptol=1e-16,gtol=1e-15,maxiter=maxit)
}else if (optim.method=="nmkb"){
    optim.arg=list(trace=NA,fnscale=NA,factr=NA,maxit=NA)
}else if (optim.method=="gensa"){
    optim.arg=list(max.time=50)
}


if ("status" %in% colnames(min.cost)){
    if (min.cost$status=="rerun_necessary"){
        foo=perform.pl(parameter.list=parm.bounds,data.matrix=all.data,pl.par.name=pl.par.name,step.size.multiplier=step.size.multiplier,optim.method=optim.method,optim.arg=optim.arg,out.file=out.file)
    }else{
        write_csv(tibble(NULL),path=out.file)
    }
}else{
    foo=perform.pl(parameter.list=parm.bounds,data.matrix=all.data,pl.par.name=pl.par.name,step.size.multiplier=step.size.multiplier,optim.method=optim.method,optim.arg=optim.arg,out.file=out.file)
}
