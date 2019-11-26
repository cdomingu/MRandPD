##This code is for Runing MR Analysis over any two Summary statistics datasets folowing the ordered format:
#RsID chr pos EA OA EAF beta se pvalue N

#RsID -> the SNP ID
#chr -> chromosome
#pos -> position 
#EA -> effect allele: allele that was selected to calculate the beta on the GWAS 
#OA -> other allele: allele that was NOT used to alculate the beta
#EAF -> allele frequency of the effect allele
#beta -> beta
#se -> standard error of the beta
#pvalue -> pvalue of the beta
#N -> sample size  

#This Script is conformed by the following steps:
#0. Arranging datasets (scripts to solve some of the following issues):
	#0.1 Dataset is lacking either RsID or chr:pos
	#0.2 Dataset has OR/TSTAT instead of Beta/SE
#1. Data Extraction and preparation
	#1.1 Clumping
	#1.2 Instrumental Variables (index SNPs) data extraction from exposure dataset
	#1.3 Instrumental Variables (index SNPs) data extraction from outcome dataset
#2. Harmonization
	#2.1 Make all the betas positive on the exposure dataset 
	#2.2 Same Effect Allele Coding Harmonization (to make sure that both the exposure and the outcome datasets have the same allele coding)
	#2.3 Calculate Rsquare and F statistic
#3. Input Generation
	#3.1 IVW
	#3.2 MRPRESSO
	#3.3 GSMR
#MR ANALYSIS
#4. IVW/MR-EGGER and MRPRESSO
#5. GSMR


###This is an example code, where the exposure is some summary statistics for intelligence (Davies,2011). 
###The user should change the word 'exposure' thoughout the code for the exposure that he or she is evaluating, i.e. TeaDrinking

#cd /data/neurogen/MRandPD/GWASsummaryStats/
#awk '{print $2,$4,$5}' snp150.txt | sed -e 's/ /:/' | sort > reference_chrposrsID
#cd /data/neurogen/MRandPD/GWASsummaryStats/Nalls2019AllSamples/
#echo "RsID chr pos EA OA EAF beta se pvalue N" > header
#join nallsEtAl2019_allSamples_allVariants_sorted /data/neurogen/MRandPD/GWASsummaryStats/reference_chrposrsID | sed -e 's/chr//g' | sed -e 's/:/ /g' | sort | awk '{print $11,$1,$2,$3,$4,$5,$6,$7,$8,$9+$10}' > temp
#cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/GWASsummaryStats/Parkinsons_GWASsummarystats.txt

#1. DATA EXTRACTION AND PREPARATION

#1.1 CLUMPING
#1.1.1 Prepare GWAS files for SNP clumping (in PLINK) by extracting the columns (SNP CHR BP P) from the GWAS Results
mkdir /data/neurogen/MRandPD/Results/PLINKClumping
mkdir /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping
cd /data/neurogen/MRandPD/GWASsummaryStats/
echo "SNP CHR BP P" > header
awk '{ if (NR!=1) print $1,$2,$3,$9}' exposure_GWASsummarystats.txt > temp
cat header temp > /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/exposure4clumping

awk '{print $1}' temp > /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/AllSNPs_exposure

# 1.1.2 Create .lsf files to run in cluster queue
#Parameters for clumping: --clump-p1 .00000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose
#this means that index snps will reach genomewide significance (.00000005) and either have a r2=0.001 with the next index SNP or be 10000 kb appart
#--clump-verbose to get the Rsquare values of the clumped SNPs with the index SNP
mkdir /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Exposure_ClumpScript/
mkdir /data/neurogen/MRandPD/Results/PLINKClumping/Exposure_Clumping/
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
cp /data/neurogen/MRandPD/Scripts/test.lsf /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Exposure_ClumpScript 
mv /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Exposure_ClumpScript/test.lsf /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Exposure_ClumpScript/exposure.clump_chr"$i".lsf
echo "module load plink/1.90b3" >> /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Exposure_ClumpScript/exposure.clump_chr"$i".lsf
echo "plink --bfile /data/neurogen/MRandPD/1000GenomsPlinkFormat/1kg_phase1_chr"$i" --extract /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/AllSNPs_exposure --clump /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/exposure4clumping --clump-p1 .00000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose --out /data/neurogen/MRandPD/Results/PLINKClumping/Exposure_Clumping/clumped.exposure.chr"$i" " >> /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Exposure_ClumpScript/exposure.clump_chr"$i".lsf
done

#1.1.3 Submit them to the queue
cd /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Exposure_ClumpScript/
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
bsub -q normal < exposure.clump_chr"$i".lsf
done 

#1.1.4 Verify that all clumping was compleated. to do so, run the following command and it whould return 22. if not, look for the ones missing and rerun them 
cd /data/neurogen/MRandPD/Results/PLINKClumping/Exposure_Clumping/
ls | grep 'log' | wc -l

#1.1.5 Extract the list of index SNPs from the clumpling results
cat clumped.exposure.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPs_exposure

#1.2 EXTRACT EXPOSURE GWAS INFO FOR CALCUATING PHENOTYPIC VARIANCE EXPALINED (Rsquare) AND F-STATISTIC
mkdir /data/neurogen/MRandPD/Results/IndexSNPs_inExposuresGWASdata/
echo "RsID chr pos EA OA EAF beta se pvalue N" > header
grep -w -f indexSNPs_exposure /data/neurogen/MRandPD/GWASsummaryStats/exposure_GWASsummarystats.txt | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/IndexSNPs_inExposuresGWASdata/exposureIVs_inExposureGWASdata

#1.3 EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET
mkdir /data/neurogen/MRandPD/Results/IndexSNPs_inParkinsonsGWASdata/
grep -w -f indexSNPs_exposure /data/neurogen/MRandPD/GWASsummaryStats/Parkinsons_GWASsummarystats.txt | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/IndexSNPs_inParkinsonsGWASdata/exposureIVs_inParkinsonsGWASdata



#2 HARMONIZATION

#2.1 MAKE ALL THE BETAS POSITIVE ON THE EXPOSURE DATASET
#It is important to convert all the betas into the same direction (positive) as it facilitates interpretation of plots and it is requiered for MR-PRESSO (see Hartwig,2016)
#Using linux shell, we'll prepare the data for beta harmonization, which will be performed on an R script. Further details on this step can be found on HarmonizationBetaPositive.R script
cd /data/neurogen/MRandPD/Results/IndexSNPs_inExposuresGWASdata/
mkdir /data/neurogen/MRandPD/Results/BetaPositiveHarmonized/
Rscript --vanilla /data/neurogen/MRandPD/Scripts/Rscripts/BetaPositiveHarmonization.R --input /data/neurogen/MRandPD/Results/IndexSNPs_inExposuresGWASdata/exposureIVs_inExposureGWASdata --exposure exposure

#2.2 SAME EFFECT ALLELE CODING HARMONIZATION
#2.2.1 After harmonizing the IV to all positive betas (BetaPositiveHarmonizarion.rmd) extract the effect allele coding from the exposure dataset. #ver que onda con R2 and F-statistic stimation
cd /data/neurogen/MRandPD/Results/BetaPositiveHarmonized/
mkdir /data/neurogen/MRandPD/Results/BetaPositiveHarmonized/ExposureGWASdata_EA/
mkdir /data/neurogen/MRandPD/Results/R2andFstatCalculation
awk '{ if (NR!=1) print $1}' /data/neurogen/MRandPD/Results/IndexSNPs_inParkinsonsGWASdata/exposureIVs_inParkinsonsGWASdata > temp
grep -w -f temp exposure_betapositive | sed -e 's/ /\t/g' | sort > /data/neurogen/MRandPD/Results/R2andFstatCalculation/exposureIVs_inExposureGWASdata
awk '{print $1,$4}' /data/neurogen/MRandPD/Results/R2andFstatCalculation/exposureIVs_inExposureGWASdata | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/BetaPositiveHarmonized/ExposureGWASdata_EA/exposureIVs_inExposureGWASdata_EA

#2.2.2 EXTRACT THE EFFECT ALLELE CODING FROM THE PD DATASET
mkdir /data/neurogen/MRandPD/Results/BetaPositiveHarmonized/ParkinsonsGWASdata_EA/
awk '{ if (NR!=1) print $1,$4}' /data/neurogen/MRandPD/Results/IndexSNPs_inParkinsonsGWASdata/exposureIVs_inParkinsonsGWASdata | sed -e 's/ /\t/g' | sort > /data/neurogen/MRandPD/Results/BetaPositiveHarmonized/ParkinsonsGWASdata_EA/exposureIVs_inParkinsonsGWASdata_EA

#2.2.3 Find which SNPs Effect Allere are different between the trait and PD datasets, for then flipping them following the same approach as before (R script)
mkdir /data/neurogen/MRandPD/Results/EAHarmonizationNallsPDdata/
mkdir /data/neurogen/MRandPD/Results/EAHarmonizationNallsPDdata/SameEAexposureAndPD/
mkdir /data/neurogen/MRandPD/Results/EAHarmonizationNallsPDdata/SNPsinPD4Flipping/
diff ExposureGWASdata_EA/exposureIVs_inExposureGWASdata_EA ParkinsonsGWASdata_EA/exposureIVs_inParkinsonsGWASdata_EA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > ../diffEA_exposurevsPD
cd /data/neurogen/MRandPD/Results/
grep -v -f diffEA_exposurevsPD IndexSNPs_inParkinsonsGWASdata/exposureIVs_inParkinsonsGWASdata > EAHarmonizationNallsPDdata/SameEAexposureAndPD/exposureIVs_inParkinsonsGWASdata_SameEA
grep -f diffEA_exposurevsPD IndexSNPs_inParkinsonsGWASdata/exposureIVs_inParkinsonsGWASdata > EAHarmonizationNallsPDdata/SNPsinPD4Flipping/exposureIVs_inParkinsonsGWASdata_EA4flip

#2.2.4 Run SameEAHarmonization.R to conduct Same Effect Allele coding harmonization
mkdir /data/neurogen/MRandPD/Results/SameEAHarmonized/
Rscript --vanilla /data/neurogen/MRandPD/Scripts/Rscripts/SameEAHarmonization.R --input /data/neurogen/MRandPD/Results/EAHarmonizationNallsPDdata/SNPsinPD4Flipping/exposureIVs_inParkinsonsGWASdata_EA4flip --exposure exposure

#2.2.5 After fliping the SNPs that needed to be fliped, join them with the ones that already had the same Effect allele coding as their corresponding exposure GWAS data.
cd /data/neurogen/MRandPD/Results/SameEAHarmonized/
mkdir /data/neurogen/MRandPD/Results/Harmonized_IndexSNPsinParkinsonsGWASdata/
cat /data/neurogen/MRandPD/Results/EAHarmonizationNallsPDdata/SameEAexposureAndPD/exposureIVs_inParkinsonsGWASdata_SameEA exposure_EAflipped | sed -e 's/ /\t/g'| sort > /data/neurogen/MRandPD/Results/SameEAHarmonized/exposureIVsinParkinsonsGWASdata_harmonized

#/data/neurogen/MRandPD/Results/SameEAHarmonized/exposureIVsinParkinsonsGWASdata_harmonized will contain Parkinson's GWAS data that will serve as the input for MR
#/data/neurogen/MRandPD/Results/R2andFstatCalculation/exposureIVs_inExposureGWASdata will contain the exposure GWAS data that will serve as the input for MR

#2.3 Calculate Rsquare and Fstatistic
Rscript --vanilla /data/neurogen/MRandPD/Scripts/Rscripts/R2andFstatistic_calculation.R --input /data/neurogen/MRandPD/Results/R2andFstatCalculation/exposureIVs_inExposureGWASdata --exposure exposure

#3 INPUT GENERATION

#INPUT GENERATION
mkdir /data/neurogen/MRandPD/Results/Input4IVW/
mkdir /data/neurogen/MRandPD/Results/Input4MRPRESSO/
mkdir /data/neurogen/MRandPD/Results/Input4GSMR/

#3.1 Input for IVW
cd /data/neurogen/MRandPD/Results/R2andFstatCalculation/
echo "bx bxse by byse snps ea oa eaf" > headerIVW
awk '{ if (NR!=1) print $1,$7,$8,$4,$5,$6}' exposureIVs_inExposureGWASdata | sort > tempIV
awk '{ if (NR!=1) print $1,$7,$8}' /data/neurogen/MRandPD/Results/SameEAHarmonized/exposureIVsinParkinsonsGWASdata_harmonized | sort > tempPD
join tempIV tempPD > temp
awk '{ print $2,$3,$7,$8,$1,$4,$5,$6}' temp > PDIV
cat headerIVW PDIV | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/Input4IVW/exposurevsPD_IVWinput
 
#3.2 Input for MRPRESSO
echo "bx bxse bxpval by byse bypval" > headerMRPRESSO
awk '{ if (NR!=1) print $1,$7,$8,$9}' exposureIVs_inExposureGWASdata  | sort > tempIV
awk '{ if (NR!=1) print $1,$7,$8,$9}' /data/neurogen/MRandPD/Results/SameEAHarmonized/exposureIVsinParkinsonsGWASdata_harmonized | sort > tempPD
join tempIV tempPD > temp
awk '{print $2,$3,$4,$5,$6,$7}' temp > PDIV
cat headerMRPRESSO 	PDIV | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/Input4MRPRESSO/exposurevsPD_MRPRESSOinput

#3.2 Input for GSMR
echo "SNP a1 a2 a1freq bx bxse bxpval bxn by byse bypval byn" > headerGSMR
awk '{ if (NR!=1) print $1,$4,$5,$6,$7,$8,$9,$10}' exposureIVs_inExposureGWASdata | sort > tempIV
awk '{ if (NR!=1) print $1,$7,$8,$9,$10}' /data/neurogen/MRandPD/Results/SameEAHarmonized/exposureIVsinParkinsonsGWASdata_harmonized  | sort > tempPD
join tempIV tempPD > IVPD
cat headerGSMR IVPD | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/Input4GSMR/exposurevsPD_GSMRinput

#4. Running IVW/MR-EGGER and MRPRESSO MR Analysis
#Run the /data/neurogen/MRandPD/Scripts/Rscripts/MRAnalysis_IVWMRPRESSO.Rmd script. In orther to use the MRPRESSO package, it requires R version below 3.5

#5. Running IVW/MR-EGGER and MRPRESSO MR Analysis
cd /data/neurogen/MRandPD/Results/Input4GSMR/
awk '{ if (NR!=1) print $1,$2}' exposurevsPD_GSMRinput > exposurevsPD_GSMRsnps.allele
module load gcta/1.91.1
gcta64 --bfile /data/neurogen/MRandPD/1000GenomsPlinkFormat/allchr --extract exposurevsPD_GSMRsnps.allele --update-ref-allele exposurevsPD_GSMRsnps.allele --recode --out gsmr_exposurevsPD
#Run the /data/neurogen/MRandPD/Scripts/Rscripts/MRAnalysis_GSMR.Rmd script.
