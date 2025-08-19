# Predicting oxygen thresholds of marine taxa to improve ecological forecasts
We used phylogenetic imputation on experimental data on taxonomically diverse marine animal taxa that related critical oxygen pressure (pO2crit) to temperature and body size. We find that pcrit declines with temperature and body size, but little similarity in the scaling across taxonomic groups.  Phylogentically imputed pcrit values led to greater predictive performance in species distribution models of groundfish in the Northeast Pacific Ocean.

# Code version 
v.2.0
# Overview of folders/files and their contents
All working R files are in subdirectory "code"
- 00-plot_taxa.R. Plots the phylogeny of species represented in the dataset and the estimation method
- 01-configure_mi_data.R.  Adds the taxonomic information by accessing the World Register of Marine Species.  Does not need to be run, as results have been generated and saved in /data/alldata_taxonomy.RDS
- 02-fit_TMB_model.R.  Runs Hierarchical model using TMB and returns many high level summaries and estimates. 
- 03-get_pcrit_for_species.R.  Code to get parameter estimates and pcrit estimates for animal taxonomic groups within and out of dataset.  Calls mane other functions including:
- - augment_taxa: adds blank taxa at all different levels for prediction of out of sample groups
- - fit_augmented_taxa: configures data and runs TMB on larger set of random effects for out of sample prediction. Set fitnew = T to rerun estimate, or F to load saved model fit
- - estimate_taxa: uses TMB model output to extract parameter estimates for taxa and estimate pcrit (returns SE of prediction).  Uses variance-covariance matrix of full precision matrix 
- 04-sdm.R.  Fits alternative species distribution models for four groundfish species, and summarizes model selection statistics using cAIC
- helper/fit_model_funs.R.  Contains numerous functions that are called in the above source files.  It is called within those source files.
- helper/plot.phylo2.funs.R.  Modified version of plot.phylo in ape package
- helper/calculate_EDF_fn.R.  Calculates empirical degrees of freedom for cAIC calculations (used by 02-fit_tmb_model.R

Files in subdirectory "analysis" are saved products as described above.  It also includes modelfit.RDS, which is a saved TMB fitted  model.

All raw data are provided in subdirectory "data".  The primary file is allmidata_aphiaxlsx.  The fields are:
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



