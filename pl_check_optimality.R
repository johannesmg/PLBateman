library(dplyr)
library(tidyr)
library(readr)
library(tibble)

args = commandArgs(trailingOnly=TRUE)


in.file=args[1]
out.file=args[2]


fit=read_csv(in.file)

fit %>%
    filter(message=="pre pl optimum") ->
    pre.pl.optimum

pre.pl.cost=pre.pl.optimum$cost[1]

fit %>%
    top_n(1,-cost) ->
    min.cost


min.cost=min.cost[1,]

spread(min.cost,parameter,parameter.value) ->
    optimal.parameters

output.status=bind_cols(tibble(status="rerun_necessary",cost.initial=pre.pl.cost,cost.final=as.numeric(min.cost[,"cost"])),optimal.parameters)

if (as.numeric(min.cost[,"cost"])>=pre.pl.cost-0.01){
    output.status$status="rerun_not_necessary"
}

print(out.file)
write_csv(output.status,path=out.file)
