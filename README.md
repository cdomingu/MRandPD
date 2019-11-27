# MRandPD
Mendelian Randomization and Parkinson's disease Project

Here, we will be using Summary Statistics from different exposure GWAS and from Parkinson's disease GWAS to determine if those exposures of interest have a causal influence over Parkinson's Disease.

We will be following the next pipeline for conducting a Mendelian Randomization study:
1.	Get the SNPs associated the exposure from its respective GWAS summary statistics.
2.	Clump them to get independent sets of index SNPs that can be used in the analysis.  (The instrumental variables)
3.	Extract those variants from the Parkinson’s Disease summary statistics.
a.	(if a variant is missing, replace it with a proxy variant)
4. Calculate the variance explained by them. Test the significance of the association between the IV and the exposure 
5.	Undergo Harmonization of Datasets
a.	Ensure that all instrumental variables are associated with the exposure on the same direction (positive), if not, flip them.
b.	Ensure that both datasets (the one from the exposure and the one from PD) are identically coded regarding the effect alleles. 
6. Apply an inverse-variance weighted (IVW) method to obtain the causal estimate.
7. Test for pleiotropy by MR-Egger regression method and check this by funnel plots.
8.	Perform MR-PRESSO as an alternative method to test for pleiotropy
9. Run the analysis under the GSMR method.
10.	Implement bi-directional MR if possible
Repeat this process for every exposure proposed.

The MRAnalysis_examplePipeline.sh is an example bash script with all the code needed to run any MR study by providing twho set of summary statistics in the specified format. The scripts BetaPositiveHarmonization.R, SameEAHarmonization.R and R2andFstatistic_calculation.R will be called thoughout the MRAnalysis_examplePipeline.sh script, in which the las step will provide the input files for the IVW, MRPRESSO and GSMR Mendelian Randomization methods.
The scripts for the MRAnalysis are intended for RStudio. The MRAnalysis_IVWMRMPRESSO.Rmd is needed to run the IVW and MRPRESSO methods. For running GSMR, first run the last lines in MRAnalysis_examplePipeline.sh and then proceed to run the script MRAnalysis_GSMR.Rmd
