treatment.model <- function(t,conc,parms){
    g1p <- parms[1]
    g2p <- parms[2]
    g1 <- parms[9]
    g2 <- parms[10]
    k01 <- parms[3]
    n1 <- parms[4]
    k02 <- parms[5]
    n2 <- parms[6]
    gg1 <- create.gamma.curve(conc,g1p,g1,k01,n1)
    gg2 <- create.gamma.curve(conc,g2p,g2,k02,n2)
    bateman.model(t,gg1,gg2,parms[7:12])
}


bateman.model <- function(t,g1,g2,parms){
    A=parms[1]
    B=parms[2]
    n=parms[5]
    t=t+parms[6]
    #make use of Kummer's transformation to prevent
                                        #machine size infinities
    g1.greater <- g1>g2
    g1.smaller <- !g1.greater
    rv <- rep(0,length(g1))
    if (sum(g1.greater)>0){
        rv[g1.greater] <- A+B*exp(-g2[g1.greater]*t[g1.greater])*t[g1.greater]^n*g1[g1.greater]^n*hyperg_1F1(a=n,b=1+n,x=-(g1[g1.greater]-g2[g1.greater])*t[g1.greater])/gamma(n+1)
    }
    if (sum(g1.smaller)>0){
        rv[g1.smaller] <- A+B*exp(-g1[g1.smaller]*t[g1.smaller])*t[g1.smaller]^n*g1[g1.smaller]^n*hyperg_1F1(a=1,b=1+n,x=(g1[g1.smaller]-g2[g1.smaller])*t[g1.smaller])/gamma(n+1)
    }
#    log2(rv)
    log2(rv+1)
}


create.gamma.curve <- function(conc,gp,g,k,n){
                                        # this should only happen for large n
                                        # and large n's should not be taken...
    if (k^n>0 & max(conc)>0){
        rv <- (gp-g)*conc^n/(k^n+conc^n)+g
    }else{
        rv <- rep(g,length(conc))
        names(rv) <- names(conc)
    }
    rv
}



load.bounds.from.file <- function(file.name){
    p0 <- read.csv(file=file.name,stringsAsFactors=FALSE)
    rownames(p0) <- p0[,"X"]
    p0 <- p0[,-1]
    list(lo.pars=as.numeric(p0["lo.pars",]),hi.pars=as.numeric(p0["hi.pars",]),start.pars=p0["start.pars",])
}

provide.measurement.and.model=function(pars,data.matrix,fn.string){
    responses <- treatment.model(data.matrix[,"t"],data.matrix[,"conc"],pars)
    w.conc=which(data.matrix[,"t"]==144 & data.matrix[,"conc"]>0 & data.matrix[,"batch"]==3)
    w.utr=which(data.matrix[,"t"]==144 & data.matrix[,"batch"]==1)
    responses[w.conc] <- responses[w.conc]-responses[w.utr]
    measured.dep <- data.matrix[,"expression"]
    measured.sd <- data.matrix[,"sd.expression"]
    list(model.exp=responses,measured.exp=measured.dep,measured.sd=measured.sd)
    }

cost.function.treated <- function(pars,data.matrix,fn.string){
    measurement.and.model=provide.measurement.and.model(pars,data.matrix,fn.string)
    sum((measurement.and.model$measured.exp-measurement.and.model$model.exp)^2/(measurement.and.model$measured.sd^2))
}

cost.function.residual.vector <- function(pars,data.matrix,fn.string){
    measurement.and.model=provide.measurement.and.model(pars,data.matrix,fn.string)
    (measurement.and.model$measured.exp-measurement.and.model$model.exp)/measurement.and.model$measured.sd
}


jacobian.cost.function <- function(pars,data.matrix){
    measured.sd <- data.matrix[,"sd.expression"]
    term.1=1/(measured.sd)
    term.2=analytical.gradient.treatment.model(data.matrix[,"t"],data.matrix[,"conc"],data.matrix[,"batch"],pars)
    term.1*term.2
}

analytical.gradient.cost.function <- function(pars,data.matrix){
    measurement.and.model=provide.measurement.and.model(pars,data.matrix,fn.string)
    term.1=2*(measurement.and.model$measured.exp-measurement.and.model$model.exp)/(measurement.and.model$measured.sd^2)
    term.2=analytical.gradient.treatment.model(data.matrix[,"t"],data.matrix[,"conc"],data.matrix[,"batch"],pars)
    as.numeric(matrix(term.1,nrow=1) %*% term.2)
}


prepare.parameters.pl=function(other.parameters,pl.parameter.value,pl.parameter.name){
    parameters=rep(0,12)
    names(parameters)=c("g1p","g2p","k01","n1","k02","n2","a","b","g1","g2","n","dt")
    pl.loc=which(names(parameters)==pl.parameter.name)
    other.loc=setdiff(1:12,pl.loc)
    parameters[pl.loc]=pl.parameter.value
    parameters[other.loc]=other.parameters
    list(parameters=parameters,other.loc=other.loc)
}

pl.interface=function(other.parameters,pl.parameter.name,pl.parameter.value,data.matrix){
    prepared.parameters=prepare.parameters.pl(other.parameters,pl.parameter.value,pl.parameter.name)
    cost.function.treated(prepared.parameters$parameters,data.matrix)
}

pl.interface.residual=function(other.parameters,pl.parameter.name,pl.parameter.value,data.matrix){
    prepared.parameters=prepare.parameters.pl(other.parameters,pl.parameter.value,pl.parameter.name)
    cost.function.residual.vector(prepared.parameters$parameters,data.matrix)
}

pl.interface.gradient=function(other.parameters,pl.parameter.name,pl.parameter.value,data.matrix){
    prepared.parameters=prepare.parameters.pl(other.parameters,pl.parameter.value,pl.parameter.name)
    -analytical.gradient.cost.function(prepared.parameters$parameters,data.matrix)[prepared.parameters$other.loc]
}

pl.interface.jacobian=function(other.parameters,pl.parameter.name,pl.parameter.value,data.matrix){
    prepared.parameters=prepare.parameters.pl(other.parameters,pl.parameter.value,pl.parameter.name)
    -jacobian.cost.function(prepared.parameters$parameters,data.matrix)[,prepared.parameters$other.loc]
}


do.treated.fitting.pl <- function(data.matrix,parameter.vector,pl.parameter.name,pl.parameter.value,optim.method,optim.arguments){
    if (optim.method=="lbfgsb"){
        yy <- optim(par=parameter.vector[["start"]],fn=pl.interface,lower=parameter.vector[["lo"]],upper=parameter.vector[["hi"]],gr=pl.interface.gradient,data.matrix=data.matrix,pl.parameter.value=pl.parameter.value,pl.parameter.name=pl.parameter.name,method="L-BFGS-B",control=optim.arguments)
    }else if (optim.method=="lm"){
        yy=nls.lm(par=parameter.vector[["start"]],fn=pl.interface.residual,lower=parameter.vector[["lo"]],upper=parameter.vector[["hi"]],jac=pl.interface.jacobian,data.matrix=data.matrix,pl.parameter.value=pl.parameter.value,pl.parameter.name=pl.parameter.name,control=optim.arguments)
    }else if (optim.method=="gensa"){
        yy <- GenSA(par=parameter.vector[["start"]],fn=pl.interface,lower=parameter.vector[["lo"]],upper=parameter.vector[["hi"]],data.matrix=data.matrix,pl.parameter.value=pl.parameter.value,pl.parameter.name=pl.parameter.name,control=optim.arguments)
    }else stop(cat("Optim method must be one of \"lbfgsb\", \"lm\", or \"gensa\""))
    yy
}


pl.wrapper.target=function(optimum.list,pl.parameter.name,extreme.val,direction,step.size,current.cost,max.cost,data.matrix,optim.method,optim.arguments){
    print(optim.arguments)
    pl.loc=which(names(optimum.list[["start"]])==pl.parameter.name)
    other.loc=setdiff(1:12,pl.loc)
    parm.bounds=lapply(optimum.list,function(x){x[other.loc]})
    par.val=log10(optimum.list[["start"]][pl.loc])
    out.vals=list()
    out.object=list()
    out.par.list=c()
    step.size=step.size*direction
    if (direction==1){
        in.bounds=par.val+1e-12<log10(extreme.val)
    }else{
        in.bounds=par.val-1e-12>log10(extreme.val)
    }
    i=1
    while (in.bounds & current.cost<max.cost){
        cat(paste(as.character(i),"...",sep=""))
        out.object[[i]]=do.treated.fitting.pl(data.matrix=data.matrix,parameter.vector=parm.bounds,pl.parameter.name=pl.parameter.name,10^par.val,optim.method,optim.arguments)
        if (optim.method=="lbfgsb"){
            out.vals[[i]]=c(out.object[[i]]$par,out.object[[i]]$value)
        }else if (optim.method=="lm"){
                out.vals[[i]]=c(out.object[[i]]$par,out.object[[i]]$deviance)
        }else if (optim.method=="nmkb"){
                out.vals[[i]]=c(out.object[[i]]$par,out.object[[i]]$value)
        }else if (optim.method=="gensa"){
                out.vals[[i]]=c(out.object[[i]]$par,out.object[[i]]$value)
        }
        parm.bounds[["start"]]=out.vals[[i]][1:11]
        out.par.list=c(out.par.list,10^par.val)
        current.cost=out.vals[[i]][12]
        par.val=par.val+step.size
        if (direction==1){
            in.bounds=par.val+1e-12<log10(extreme.val)
        }else{
            in.bounds=par.val-1e-12>log10(extreme.val)
        }
        i=i+1
    }
    if (length(out.par.list)>0){
        messages=as.character(lapply(out.object,"[[","message"))
        out.vv=bind_cols(tibble(g1=out.par.list),as_tibble(do.call("rbind",out.vals)))
        out.vv$message=messages
        colnames(out.vv)=c("parameter.value",names(optimum.list[["start"]])[other.loc],"cost","message")
    }else{
        out.vv=NULL
        }
    out.vv
}


provide.starting.parameters <- function(gene.name,seed.no,no.reps,no.pars,par.names,bounds.frame,random.gen="randomLHS"){
    fran <- get(random.gen)
    set.seed(seed.no)
    random.pars <- fran(no.reps,no.pars)
    random.pars <- lapply(1:no.pars,function(y){10^(log10(bounds.frame$lo.observed[y])+(log10(bounds.frame$hi.observed[y])-log10(bounds.frame$lo.observed)[y])*random.pars[,y])})
    random.pars <- as_tibble(do.call("cbind",random.pars))
    colnames(random.pars) <- par.names
    bind_cols(tibble(gene=rep(gene.name,nrow(random.pars))),random.pars)
}    



perform.pl=function(parameter.list,data.matrix,pl.par.name,step.size.multiplier,optim.method,optim.arg,out.file){

    yy <- optim(par=parameter.list[["start"]],fn=cost.function.treated,lower=parameter.list[["lo"]],upper=parameter.list[["hi"]],gr=function(x,...){-analytical.gradient.cost.function(x,...)},data.matrix=data.matrix,method="L-BFGS-B",control=list(maxit=20000,fnscale=1e6,factr=1e2))

    ix=which(par.names==pl.par.name)

    parameter.list$start=yy$par
    names(parameter.list$start)=par.names

    par.lower.bound=parm.bounds$lo[ix]
    par.higher.bound=parm.bounds$hi[ix]

    cost.increase=qchisq(0.95,1)
    step.size=(log10(par.higher.bound)-log10(par.lower.bound))/step.size.multiplier

    enough.points=FALSE
    while (!enough.points){
        all.pl.r=pl.wrapper.target(parameter.list,pl.par.name,extreme.val=par.higher.bound,direction=1,step.size=step.size,current.cost=yy$value,max.cost=yy$value+cost.increase,data.matrix=data.matrix,optim.method,optim.arguments = optim.arg)
        print(length(all.pl.r))
        if (is.null(all.pl.r)){
            enough.points=TRUE
        }else{        
            if (nrow(all.pl.r)>=10){
                enough.points=TRUE
            }else{
                step.size=step.size/5
            }
        }
    }

    enough.points=FALSE
    step.size=(log10(par.higher.bound)-log10(par.lower.bound))/step.size.multiplier

    while (!enough.points){
        all.pl.l=pl.wrapper.target(parameter.list,pl.par.name,extreme.val=par.lower.bound,direction=-1,step.size=step.size,current.cost=yy$value,max.cost=yy$value+cost.increase,data.matrix=data.matrix,optim.method,optim.arguments = optim.arg)
        if (is.null(all.pl.l)){
            enough.points=TRUE
        }else{        
            if (nrow(all.pl.l)>=10){
                enough.points=TRUE
            }else{
                step.size=step.size/5
            }
        }
    }

    out.data=bind_rows(all.pl.l,all.pl.r)
    out.data$method=optim.method
    out.data$gene=current.gene
    out.data$parameter=pl.par.name
    
    par.frame=enframe(yy$par[-ix]) %>%
        spread(name,value) %>%
        mutate(parameter=pl.par.name,parameter.value=yy$par[ix],cost=yy$value,gene=current.gene,method=optim.method,message="pre pl optimum")
    out.data=bind_rows(par.frame,out.data)
    
    if (optim.method=="lbfgsb"){
        out.data=bind_cols(out.data,as_tibble(optim.arg)[rep(1,nrow(out.data)),])
    }else if (optim.method=="lm"){
        out.data=bind_cols(out.data,as_tibble(optim.arg[c("ftol","ptol","gtol","maxiter")])[rep(1,nrow(out.data)),])
    }else if (optim.method=="nmkb"){
        out.data=bind_cols(out.data,as_tibble(optim.arg)[rep(1,nrow(out.data)),])
    }else if (optim.method=="gensa"){
        out.data=bind_cols(out.data,as_tibble(optim.arg)[rep(1,nrow(out.data)),])
    }
    write_csv(out.data,path=out.file)
}


provide.data=function(data.file,current.gene){
    all.data=read_csv(data.file,guess_max=3e4)
    all.data=all.data %>%
        ungroup %>%
        filter(gene == current.gene)

    all.data <- all.data[,c("t","mny","sdy","conc","batch","gene")]
    colnames(all.data) <- c("t","expression","sd.expression","conc","batch","gene")

    all.data %>%
        select(-gene)->
        all.data

    all.data=as.matrix(all.data)
}
