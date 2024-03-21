# Remarkable similarity of oxygen tolerance across marine taxa when standardized for temperature and body size
# Summary of the study
We used phylogenetic imputation on experimental data on taxonomically diverse marine taxa that related critical oxygen pressure (pO2crit) to temperature and body size.  We found that once we accounted for differences in body size, temperatures, and the source paper, the data point to very similar oxygen partial pressures (generally between 3 - 5 kPa) at which oxygen demand exceeds supply.  

# Code version 
v.1.0
# Overview of folders/files and their contents
All working R files are in subdirectory "code"

- 01-configure-mi-data.R.  Adds the taxonomic information by accessing the World Register of Marine Species.  Does not need to be run, as results have been generated and saved in /data/alldata_taxonomy.RDS
- 02-fit-individual-species.R. Runs linear models individually for each species.  Does not need to be run, as results have been generated and saved in analysis/species_estimates.RDS
- 03-fit_TMB_model.R.  Runs Hierarchical model using TMB, and plots results for taxonomic groups. Also includes summary of data coverage by taxonomic group.
- 04-fit_TMB_stan.R.  Runs the hierarchical model using tmbstan.  Does not need to be run, as results have been generated and saved in analysis/mcmcoutput.RDS
- 05-sim_new_taxa.R.  Take the output of tmbstan, and then simulates a random species within each specified taxonomic level
- 06-summarize_sims.R. Code used to summarize results, organized by taxonomic grouping and used to create supplemental table
- fit_model_funs.R.  Contains numerous functions that are called in the above source files.  It is called within those source files.

Files in subdirectory "analysis" are saved products as described above.  It also includes modelfit.RDS, which is a saved TMB fitted  model, and modelfit_stan.RDS, which is another saved TMB model where the model code is modified to permit successful MCMC integration.

All raw data are provided in subdirectory "data".  The primary file is allmidata_csv.  The fields are:
- scientific.name: Genus species, updated January 2023; https://www.marinespecies.org/index.php
- Temp: experimental temperature (degrees C)
- W: individual (or mean individual) body size (grams)
- Pcrit: estimated critical oxygen partial pressure (kPa)
- AphiaID: Aphia ID for species, genus, family, or order: https://www.marinespecies.org/index.php
- lowest.taxon: lowest taxonomic level for taxonomic group.  Usually equals "species"
- Source: Short verson of paper citation or data source

# Instructions for users to run the software (e.g. explain the project workflow and any configuration parameters of your software)
See workflow instructions above.  To extract simulations for a taxa, go to the lookup_taxa sub directory and open file "lookup_taxa.R".  There you can specify the taxonomic group, and source the file to see simulation medians, standard deviations, and variance - covariance matrix of metabolic index tratis


# Links to protocols.io or equivalent methods repositories, where applicable
NA



