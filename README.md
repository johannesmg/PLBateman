These R scripts compute the profile likelihood for fit parameters for the article "Kinetic modeling of stem cell transcriptome dynamics to identify regulatory modules of normal and disturbed neuroectodermal differentiation". Here, the GPL (>=2) applies.
 
The Snakefile describes the workflow and can be executed using snakemake, downloadable from https://github.com/snakemake/snakemake.

Before the Snakefile can be executed, fits need to be computed using the code in the repository https://github.com/johannesmg/BatemanDiff

Additionally, the following files are required:

scaled_expression_data_file.csv: This file can be created using the code in the repository https://github.com/johannesmg/BatemanDiff
parameter_boundary_file.csv: A file with suggestions for the parameter boundaries can be found in the repository https://github.com/johannesmg/BatemanDiff and should be identical to the one used for computing the fits

The configuration of the local directories has to be entered into the config.yml file.