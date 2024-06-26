# Yeast_mutation_effects_across_strains_and_environments
This repository contains all relevant code for completing analyses and producing figures for the paper 'Environment-independent distribution of mutational effects emerges from microscopic epistasis'. To use this code, you will need an updated version of R. *NOTE: All Files and Scripts here are preliminary. Finalized data and scripts will be uploaded upon the manuscript's acceptance*. To properly run this code ensure:
1. You have installed all packages listed at the top of the script.
2. You update the working directory (and, optionally, the path to save files) on lines 38-39. The code will only run if your working directory has all relevant files.

Unfortunately, the confirmation experiment BC count file is too large to upload here. It is available upon request and I have kept the code to visualize it in the script but commented it out so as not to disrupt a full run of the other figures. 

The file “CoreFigureCode.R”, when run with the relevant data files in the working directory, will produce all main text and supplemental figures. Note that the following files are required and readily available in this repository):
-	**df_s_ests_wGR** : Contains the average fitness effect estimates for all mutations in all strains and environments where we had sufficient measurements. Column names are as follows
-	  Strain: Strain names, consistent with Johnson et al 2019
-	  Env: environment name, as ‘Temp’ ‘SC’ ‘pH’
-	  Mut_ID: Numerical identifier corresponding to each mutation. Gene names for each can be found in Mut_fullinfo_use data frame
-	  avgS: Average fitness effect of each mutation
-	  varS: variance in the fitness effect
-	  seS: standard error of the fitness effect
-	  Mean_GR: Estimated growth rate for the corresponding strain in the corresponding environment
-	  Std_err: Standard error of the growth rate estimate
-	  mutCall: call of each mutation’s effect
-	  minCI: bottom of the 99% confidence interval range for the mutation’s effect
-	  maxCI: top of the 99% confidence interval range for the mutation’s effect
-	**KRE33_pres_abs**: Contains indicators of the allelic variant (0 or 1; BY vs RM) at each focal QTL locus (indicated as column names) 
-	**Mut_fullInfo_use**: Maps mutation numerical IDs to the gene they are in or nearby
-	**JohnsonDFEmeanDat** – can be created given data from Johnson et al Science 2019: simply need the strain names and the corresponding parameters: DFEmean, DFEvar, DFEskew , DFEmean_min, DFEmean_max, DFEvar_min, DFEvar_max, DFEskew_min, DFEskew_max. I have pre-processed this freely available data and uploaded the file here for ease of use.
-	  However: The growth rates provided in Johnson et al 2019 are relative growth rate metrics. We convert these to exponential growth rate using the relationship between relative and exponential growth given in their supplemental figure S9 and detailed in the supplementary materials of the manuscript accompanying this code.
-	**allBCMatches**: Contains raw data of barcode counts for the confirmation experiment. *Available upon request - smardell@ucsd.edu*
-	**transferODtable**: Contains information to convert data file names to relevant sample information for the confirmation experiment


Note that many sections depend on data files created in other sections, so it is best run all together, Total runtime is less than 10 minutes (on Mac 16GB Ram M1 chip laptop). 
