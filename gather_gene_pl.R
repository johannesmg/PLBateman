library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library (optparse)

par.names <- c("g1p","g2p","k01","n1","k02","n2","a","b","g1","g2","n","dt")


option_list <- list ( make_option (c("-f","--filelist"),default="blah.txt", 
                                   help="comma separated list of files (default %default)")
                     )

parser <-OptionParser(option_list=option_list)
arguments <- parse_args (parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args

myfilelist <- strsplit(opt$filelist, ",")

print(myfilelist)
print(args)

out.file=args[1]

fit=bind_rows(lapply(unlist(myfilelist),read_csv,col_names=T))

write_csv(fit,path=out.file)
