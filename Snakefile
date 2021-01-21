import pandas as pd
from collections import Counter

configfile: 'config.yml'

EXPRESSIONDATA=config['expressiondata']
OUTTREATED=config['outtreated']
PARAMETERFILE=config['parameterfile']
GENSATIME=config['gensatime']
NOCORESTREATEDFIT=config['nocorestreatedfit']
OUTPL=config['outpl']
OUTSECOND=config['outsecond']
INPUTPL=config['inputpl']
OPTIMMETHOD=config['optimmethod']
MULTIPLIER=config['multiplier']
FNSCALE=config['fnscale']
FACTR=config['factr']
MAXIT=config['maxit']
PARAMETER=["a","b","g1","g2","k01","n1","k02","n2","g1p","g2p","n","dt"]

INPUTMETHOD=["lbfgsb","lm","gensa","nmkb"]

ll=len(INPUTMETHOD)

(genes,mtd,fnsc,fct) = glob_wildcards(OUTTREATED+"/final_2_out_{gene}_{method}_{fnscale}_{factr}.csv")
cts=Counter(genes)
genes=list({x: count for x, count in cts.items() if count >= ll})

print(len(genes))

rule all:
  input: [expand("{outsecond}/pl_outgene_{gene}.csv",gene=genes,outsecond=OUTSECOND)]

rule interface:
  input: F=expand('{outtreated}/final_2_out_{{gene}}_{method}_1e6_{factr}.csv',outtreated=OUTTREATED,method=INPUTMETHOD,factr=FACTR)
  output: INPUTPL+'/pl_in_{gene}.csv'
  threads: 1
  resources: 
  run:
      file_lst = ",".join(map(str, input.F))
      shell('Rscript interface_to_fit_result.R -f {file_lst} {output}')

rule profile_likelihood:
  input: [EXPRESSIONDATA,INPUTPL+'/pl_in_{gene}.csv']
  output: OUTPL+'/pl_out_{parameter}_{gene}_{method}_{factr}.csv'
  threads: 1
  resources: 
    mem_per_cpu=2000
  message: "--- Fitting gene treated."
  shell: 'Rscript perform_pl.R {input[0]} {input[1]} {output} {PARAMETERFILE} {wildcards.parameter} {MULTIPLIER} {threads} {wildcards.method} {FNSCALE} {wildcards.factr} {MAXIT}'

rule gather_gene:
  input: F=expand('{outpl}/pl_out_{parameter}_{{gene}}_{method}_{factr}.csv',outpl=OUTPL,method=OPTIMMETHOD,factr=FACTR,parameter=PARAMETER)
  output: OUTPL+'/pl_outgene_{gene}.csv'
  threads: 1
  message: "--- Checking cost."
  run:
      file_lst = ",".join(map(str, input.F))
      shell('Rscript gather_gene_pl.R -f {file_lst} {output}')
      shell('rm {input.F}')

rule check_optimality:
  input: OUTPL+'/pl_outgene_{gene}.csv'
  output: OUTPL+'/checked_optimality_{gene}.csv'
  threads: 1
  message: "--- Checking cost."
  shell: 'Rscript pl_check_optimality.R {input} {output}'

rule profile_likelihood_second:
  input: [EXPRESSIONDATA,OUTPL+'/checked_optimality_{gene}.csv']
  output: OUTSECOND+'/pl_out_{parameter}_{gene}_{method}_{factr}.csv'
  threads: 1
  resources: 
    mem_per_cpu=2000
  message: "--- Fitting gene treated."
  shell: 'Rscript perform_pl.R {input[0]} {input[1]} {output} {PARAMETERFILE} {wildcards.parameter} {MULTIPLIER} {threads} {wildcards.method} {FNSCALE} {wildcards.factr} {MAXIT}'

rule gather_gene_final:
  input: F=expand('{outsecond}/pl_out_{parameter}_{{gene}}_{method}_{factr}.csv',outsecond=OUTSECOND,method=OPTIMMETHOD,factr=FACTR,parameter=PARAMETER)
  output: OUTSECOND+'/pl_outgene_{gene}.csv'
  threads: 1
  message: "--- Checking cost."
  run:
      file_lst = ",".join(map(str, input.F))
      shell('Rscript gather_gene_pl.R -f {file_lst} {output}')
      shell('rm {input.F}')
