# Predicting oxygen tolerance of marine fauna through a taxonomic hierarchical analysis
# Summary of the study
We used phylogenetic imputation on experimental data on taxonomically diverse marine taxa that related critical oxygen pressure (pO2crit) to temperature and body size. We find that pcrit declines with temperature and body size, but little similarity in the scaling across taxonomic groups.  Phylogentically imputed pcrit values led to greater predictive performance in species distribution models of groundfish in the Northeast Pacific Ocean.

# Code version 
v.1.1
# Overview of folders/files and their contents
All working R files are in subdirectory "code"
- 00-plot_taxa.R. Plots the phylogeny of species represented in the dataset and the estimation method
- 01-configure_mi_data.R.  Adds the taxonomic information by accessing the World Register of Marine Species.  Does not need to be run, as results have been generated and saved in /data/alldata_taxonomy.RDS
- 02-fit_alternative_models.R.  Fits different models with different top taxonomic levels; Class, Order, Family, and Genus.  Prediction performance is maximized when setting Family as the top level.
- 03-fit_TMB_model.R.  Runs Hierarchical model using TMB, using Class as the topmost taxonomic level, and plots group-level mean trait values. 
- 04-fit_TMB_stan.R.  Runs the hierarchical model using tmbstan.  Does not need to be run, as results have been generated and saved in analysis/mcmcoutput.RDS
- 05-sim_new_taxa.R.  Take the output of tmbstan, and then simulates a random species within each specified taxonomic level
- 06-summarize_sims.R. Code used to summarize results, organized by taxonomic grouping and used to create supplemental table
- 07-sdm.R.  Fits alternative species distribution models for four groundfish species, and summarizes model selection statistics using AIC
- helper/fit_model_funs.R.  Contains numerous functions that are called in the above source files.  It is called within those source files.
- helper/plot.phylo2.funs.R.  Modified version of plot.phylo in ape package
- helper/calculate_EDF_fn.R.  Calculates empirical degrees of freedom for cAIC calculations (uased by 02-fit_alternative_models.R

Files in subdirectory "analysis" are saved products as described above.  It also includes modelfit.RDS, which is a saved TMB fitted  model, and modelfit_stan.RDS, which is another saved TMB model where the model code is modified to permit successful MCMC integration.

All raw data are provided in subdirectory "data".  The primary file is allmidata_xlsx.  The fields are:
- scientific.name: Genus species, updated June 2025; https://www.marinespecies.org/index.php
- Temp: experimental temperature (degrees C)
- W: individual (or mean individual) body size (grams)
- Pcrit: estimated critical oxygen partial pressure (kPa)
- AphiaID: Aphia ID for species, genus, family, or order: https://www.marinespecies.org/index.php
- lowest.taxon: lowest taxonomic level for taxonomic group.  Usually equals "species"
- Source: Paper citation

# Instructions for users to run the software (e.g. explain the project workflow and any configuration parameters of your software)
See workflow instructions above.


# Links to protocols.io or equivalent methods repositories, where applicable
NA



