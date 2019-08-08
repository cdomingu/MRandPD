# MRandPD
Mendelian Randomization and Parkinson's disease Project

Here, we will be using Summary Statistics from different exposure GWAS and from Parkinson's disease GWAS to determine if those exposures of interest have a causal influence over Parkinson's Disease.

The exposures that have been included in the analysis thus far are:
1. Smoking (age of initiation, cigarrettes per day, ever vs never smoking, former vs current smoker).
2. Sleeping (excess daytime sleepiness, insomnia, sleep duration)
3. Drinking (drinks per week)
4. Educational Attainment (years of education)

We will be following the next pipeline for conducting a Mendelian Randomization study:
1.	Get the SNPs associated the exposure from its respective GWAS summary statistics.
2.	Clump them to get independent sets of index SNPs that can be used in the analysis.  (The instrumental variables)
3.	Calculate the variance explained by them. 
4.	Test the significance of the association between the instrument and the exposure using F statistic.
5.	Extract those variants from the Parkinson’s Disease summary statistics.
a.	(if a variant is missing, replace it with a proxy variant)
6.	Undergo Harmonization of Datasets
a.	Ensure that all instrumental variables are associated with the exposure on the same direction (positive), if not, flip them.
b.	Ensure that both datasets (the one from the exposure and the one from PD) are identically coded regarding the effect alleles. 
7.	Calculate the Wald estimator for each of the individual SNPs by dividing the log-OR for PD and dividing it by the mean difference of the SB for each SNP (obtained by dividing the variant-exposure association over the variant-outcome association). (log-OR will be obtained through PRS). 
8.	Apply an inverse-variance weighted (IVW) method to obtaining the linear regression. 
9.	Test for heterogeneity by MR-Egger regression method to test for bias due to horizontal pleiotropy, and check this by funnel plots. 
10.	Perform MR-PRESSO as an alternative method to test for pleiotropy
11.	Do power calculations.
12.	Implement bi-directional MR if possible
Repeat this process for every exposure proposed.


