# metabolic_index
Code for running phylogenetic trait analysis on metabolic index parameters.

All R files to run are in subdirectory "code"
- 01-configure-mi-data.R.  Adds the taxonomic information by accessing the World Register of Marine Species.  Does not need to be run, as results have been generated and saved in /data/alldata_taxonomy.RDS
- 02-fit-individual-species.R. Runs linear models individually for each species.  Does not need to be run, as results have been generated and saved in analysis/species_estimates.RDS
- 03-fit_TMB_model.R.  Runs Hierarchical model using TMB, and plots results for taxonomic groups. Also includes summary of data coverage by taxonomic group.
- 04-fit_TMB_stan.R.  Runs the hierarchical model using tmbstan.  Does not need to be run, as results have been generated and saved in analysis/mcmcoutput.RDS
- fit_model_funs.R.  Contains numerous functions that are called in the above source files.  It is called within those source files.

Files in subdirectory "analysis" are saved products as described above.  It also includes modelfit.RDS, which is a saved TMB fitted  model, and modelfit_stan.RDS, which is another saved TMB model where the model code is modified to permit successful MCMC integration.
