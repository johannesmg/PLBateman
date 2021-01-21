library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library (optparse)

option_list <- list ( make_option (c("-f","--filelist"),default="foo.txt", 
                                   help="comma separated list of files (default %default)")
                     )

parser <-OptionParser(option_list=option_list)
arguments <- parse_args (parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args

myfilelist <- strsplit(opt$filelist, ",")

out.file=args[1]

fit=bind_rows(lapply(unlist(myfilelist),read_csv,col_names=T))

fit %>%
    top_n(1,-cost) ->
    min.cost

write_csv(min.cost,path=out.file)
