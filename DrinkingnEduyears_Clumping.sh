##This code is for obtainig the index SNPs (IV) by clumping from DRINKING and EDUCATION YEARS GWAS data
##We need to clump SNPs to ensure that all the index SNPs are independent frome one another
##After that, we will extract those index SNPs data from the PD GWAS, an calculate the Rsquare an F-statistic from the instruments
##We are using the summary statistics for drinking from Liu,2019 downloaded from https://conservancy.umn.edu/handle/11299/201564
##We are using the summary statistics for eduyears from Okbay,2016 downloaded from https://www.thessgac.org/data

##DRINKING

#2. PREPARE GWAS FILES FOR SNP CLUMPING (in PLINK)
# EXTRACTION OF COLUMNS (SNP EA NEA P BETA) from the GWAS results
cd /data/neurogen/MRandPD/GWASsummaryStats/drinking/
echo "SNP CHR BP P" | sed -e 's/ /\t/g'  > header
awk '{ if (NR!=1) print $3,$1,$2,$8}' DrinksPerWeek.txt | sed -e 's/ /\t/g' > temp
cat header temp > /data/neurogen/MRandPD/Results/Clumping/Data4Clump/drinksperweek

# 3.1 COLLATE ALL SNPS FROM GWAS
cd /data/neurogen/MRandPD/Results/Clumping/Data4Clump/
awk '{ if (NR != 1) print $1}' drinksperweek >> allSnpRsGWASdrinking
# 3.3 SORT RS NUMBERS AND REMOVE DUPLICATED SNPS
sort allSnpRsGWASdrinking | uniq > allSnpRsUniquedrinking
# 3.4 GET RID OF INDELS (THEY DO NOT HAVE AN RS NUMBER)
grep rs allSnpRsUniquedrinking > allSnpRsUniqueNoIndelsdrinking

# Using the 1000Genomes data in Plink format downloaded from: http://www.cog-genomics.org/plink/1.9/resources

# 4. CREATE .lsf FILES TO RUN IN CLUSTER QUEUE
#Parameters for clumping: --clump-p1 .00000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose
#this mean that index snps will reach genomewide significance (.00000005) and either have a r2=0.001 with the next index SNP or be 10000 kb appart
#--clump-verbose to ger the Rsquare values of the clumped SNPs with the index SNP
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
cp /data/neurogen/MRandPD/Scripts/test.lsf /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/DrinkingClumpScrptpvalmin8Verbose 
mv /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/DrinkingClumpScrptpvalmin8Verbose/test.lsf /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/DrinkingClumpScrptpvalmin8Verbose/drinksperweek.clump_chr"$j".lsf
echo "module load plink/1.90b3" >> /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/DrinkingClumpScrptpvalmin8Verbose/drinksperweek.clump_chr"$j".lsf
echo "plink --bfile /data/neurogen/MRandPD/1000GenomsPlinkFormat/1kg_phase1_chr"$j" --extract /data/neurogen/MRandPD/Results/Clumping/Data4Clump/allSnpRsUniqueNoIndelsdrinking --clump /data/neurogen/MRandPD/Results/Clumping/Data4Clump/drinksperweek --clump-p1 .00000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose --out /data/neurogen/MRandPD/Results/Clumping/DrinkingClumppvalmin8Verbose/clumped.drinksperweek.chr"$j" " >> /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/DrinkingClumpScrptpvalmin8Verbose/drinksperweek.clump_chr"$j".lsf
done

# SUBMIT THEM TO THE QUEUE
cd /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/DrinkingClumpScrptpvalmin8Verbose
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
bsub -q normal < drinksperweek.clump_chr"$j".lsf
done 

#5.EXTRACT THE LIST OF INDEXSNPS FROM CLUMPING RESULTS
cd /data/neurogen/MRandPD/Results/Clumping/DrinkingClumppvalmin8Verbose/
cat clumped.drinksperweek.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPSdrinksperweekV

#6.EXTRACT EXPOSURE GWAS INFO (specially BETA, MAF,SE) FOR CALCUATING PHENOTYPIC VARIANCE EXPALINED (Rsquare) AND F-STATISTIC
echo "CHROM POS RSID REF ALT AF STAT PVALUE BETA SE N EFFECTIVE_N Number_of_Studies ANNO ANNOFULL" | sed -e 's/ /\t/g'  > header
grep -w -f indexSNPSdrinksperweekV /data/neurogen/MRandPD/GWASsummaryStats/drinking/DrinksPerWeek.txt | sort -n -k 2 | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticprePD/drinksperweekforR2calc

#7. EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET
echo "MarkerName Allele1 Allele2 Effect StdErr P-value NStudies HetISq" | sed -e 's/ /\t/g'  > header
grep -w -f indexSNPSdrinksperweekV /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/drinksperweekSNPinPDdataV

#8.GET MISSING VARIANTS (INDEX SNPS THAT WEREN'T FOUND ON PD GWAS DATASET)
cd /data/neurogen/MRandPD/Results/HarmonizationVerbose/
awk '{ if (NR!=1) print $1}' drinksperweekSNPinPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/Clumping/DrinkingClumppvalmin8Verbose/indexSNPSdrinksperweekV > missingindexSNPdrinksperweekV

#9.IDENTIFY PROXY SNPS FROM MISSING VARIANT (THE ONES WITH Rsquare > 0.9 WITH INDEX SNP)
cd /data/neurogen/MRandPD/Results/Clumping/DrinkingClumppvalmin8Verbose/
cat clumped.drinksperweek.chr*.clumped > temp
sed -e 's/     / /g' temp| sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > prueba
awk -F'[ ]' '$4>0.9 || $5==1' prueba | grep 'rs' | grep -v "^\s[0-9]" > rsqSNPsdrinksperweek

#10."MANUALLY" DETERMINE HOW MANY PROXY SNPS WERE FOUND DURING THE CLUMPING FOR THE MISSING SNPS
grep -w -A10 'everymissingsnip' rsqSNPsdrinksperweek

#11.EXTRACT THE PROXY SNPS THAT WERE FOUND FOR ANY OF THE MISSING SNPs
grep -A27 'rs34704785' rsqSNPsdrinksperweek | grep -v INDEX | cut -d ' ' -f2,4,7 > temp
awk '$2 > 0.95 && $3 < .000001' temp | cut -d ' ' -f1 > rs34704785RsqSNPdrinksperweek
grep -A178 'rs113589236' rsqSNPsdrinksperweek | grep -v INDEX | cut -d ' ' -f2,4,7 > temp
awk '$2 > 0.97 && $3 < .00000000000000002' temp | cut -d ' ' -f1 > rs113589236RsqSNPdrinksperweek
grep -A3 'rs76217384' rsqSNPsdrinksperweek | grep -v INDEX | cut -d ' ' -f2 > rs76217384RsqSNPdrinksperweek
grep -A1 'rs113909752' rsqSNPsdrinksperweek | grep -v INDEX | cut -d ' ' -f2 > rs113909752RsqSNPdrinksperweek

#12.EXTRAXT THIS NEW PROXY SNPs FROM THE PD GWAS DATASET, IF THERE ARE FOUND
grep -w -f rs34704785RsqSNPdrinksperweek /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs34704785RsqSNPdrinksperweekPD
grep -w -f rs113589236RsqSNPdrinksperweek /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs113589236RsqSNPdrinksperweekPD
grep -w -f rs76217384RsqSNPdrinksperweek /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs76217384RsqSNPdrinksperweekPD
grep -w -f rs113909752RsqSNPdrinksperweek /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs113909752RsqSNPdrinksperweekPD

#13.JOIN ALL INDEX SNPS FOUND ON THE PD GWAS DATASET
cat drinksperweekSNPinPDdataV rs34704785RsqSNPdrinksperweekPD rs113589236RsqSNPdrinksperweekPD rs76217384RsqSNPdrinksperweekPD rs113909752RsqSNPdrinksperweekPD > fulldrinksperweekSNPinPDdataV

#14.EXTRACT THE INFO FOR CALCULATING Rsquare AND F-STATISTIC FROM INDEX SNPS ALSO PRESENT IN THE PD GWAS DATASET 
echo "CHROM POS RSID REF ALT AF STAT PVALUE BETA SE N EFFECTIVE_N Number_of_Studies ANNO ANNOFULL" | sed -e 's/ /\t/g'  > header
awk '{ if (NR!=1) print $1}' fulldrinksperweekSNPinPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/drinking/DrinksPerWeek.txt  | awk '$8<.000001' | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/drinksperweekforR2calcV


##EDUCATIONAL ATTAINMENT

#2. PREPARE GWAS FILES FOR SNP CLUMPING (in PLINK)
# EXTRACTION OF COLUMNS (SNP EA NEA P BETA) from the GWAS results
# eduyears
echo "SNP CHR BP P" | sed -e 's/ /\t/g'  > header
awk '{ if (NR!=1) print $1,$2,$3,$9}' EduYears_Main.txt | sed -e 's/ /\t/g' > temp
cat header temp > /data/neurogen/MRandPD/Results/Clumping/Data4Clump/eduyears

# 3.1 COLLATE ALL SNPS FROM GWAS
cd /data/neurogen/MRandPD/Results/Clumping/Data4Clump/
awk '{ if (NR != 1) print $1}' eduyears >> allSnpRsGWASeduyears
# 3.3 SORT RS NUMBERS AND REMOVE DUPLICATED SNPS
sort allSnpRsGWASeduyears | uniq > allSnpRsUniqueeduyears
# 3.4 GET RID OF INDELS (THEY DO NOT HAVE AN RS NUMBER)
grep rs allSnpRsUniqueeduyears > allSnpRsUniqueNoIndelseduyears

# Using the 1000Genomes data in Plink format downloaded from: http://www.cog-genomics.org/plink/1.9/resources

#4. CREATE .lsf FILES TO RUN IN CLUSTER QUEUE
#Parameters for clumping: --clump-p1 .00000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose
#for identifying both genomewide significant that are independent form one another
#this mean that intex snps will reach significance (.00000005) and either have a r2=0.001 with the next index SNP or be 10000 kb appart
#--clump-verbose to ger the Rsquare values of the clumped SNPs with the index SNP
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
cp /data/neurogen/MRandPD/Scripts/test.lsf /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/EduyearsClumpScrptpvalmin8Verbose 
mv /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/EduyearsClumpScrptpvalmin8Verbose/test.lsf /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/EduyearsClumpScrptpvalmin8/eduyears.clump_chr"$j".lsf
echo "module load plink/1.90b3" >> /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/EduyearsClumpScrptpvalmin8Verbose/eduyears.clump_chr"$j".lsf
echo "plink --bfile /data/neurogen/MRandPD/1000GenomsPlinkFormat/1kg_phase1_chr"$j" --extract /data/neurogen/MRandPD/Results/Clumping/Data4Clump/allSnpRsUniqueNoIndelseduyears --clump /data/neurogen/MRandPD/Results/Clumping/Data4Clump/eduyears --clump-p1 .00000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose --out /data/neurogen/MRandPD/Results/Clumping/EduyearsClumppvalmin8Verbose/clumped.eduyears.chr"$j" " >> /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/EduyearsClumpScrptpvalmin8Verbose/eduyears.clump_chr"$j".lsf
done

# SUBMIT THEM TO THE QUEUE
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
bsub -q normal < eduyears.clump_chr"$j".lsf
done 

#5.EXTRACT THE LIST OF INDEXSNPS FROM CLUMPING RESULTS
cd /data/neurogen/MRandPD/Results/Clumping/EduyearsClumppvalmin8Verbose/
cat clumped.eduyears.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPSeduyearsV

#6.EXTRACT EXPOSURE GWAS INFO (specially BETA, MAF,SE) FOR CALCUATING PHENOTYPIC VARIANCE EXPALINED (Rsquare) AND F-STATISTIC
echo "MarkerName CHR POS A1 A2 EAF Beta SE Pval" | sed -e 's/ /\t/g'  > header
grep -w -f indexSNPSeduyearsV /data/neurogen/MRandPD/GWASsummaryStats/education/EduYears_Main.txt | sort -n -k 2 | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticprePD/eduyearsforR2calc

#7. EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET
cd /data/neurogen/MRandPD/Results/Clumping/EduyearsClumppvalmin8Verbose/
echo "MarkerName Allele1 Allele2 Effect StdErr P-value NStudies HetISq" | sed -e 's/ /\t/g'  > header
grep -w -f indexSNPSeduyearsV /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/eduyearsSNPinPDdataV

#8. GET MISSING VARIANTS (INDEX SNPS THAT WEREN'T FOUND ON PD GWAS DATASET)
cd /data/neurogen/MRandPD/Results/HarmonizationVerbose/
awk '{ if (NR!=1) print $1}' eduyearsSNPinPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/Clumping/EduyearsClumppvalmin8Verbose/indexSNPSeduyearsV > missingindexSNPeduyearsV

#9.IDENTIFY PROXY SNPS FROM MISSING VARIANT (THE ONES WITH Rsquare > 0.9 WITH INDEX SNP)
cd /data/neurogen/MRandPD/Results/Clumping/EduyearsClumppvalmin8Verbose/
cat clumped.eduyears.chr*.clumped > temp
sed -e 's/     / /g' temp| sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > prueba
awk -F'[ ]' '$4>0.9 || $5==1' prueba | grep 'rs' | grep -v "^\s[0-9]" > rsqSNPseduyears

#10."MANUALLY" DETERMINE HOW MANY PROXY SNPS WERE FOUND DURING THE CLUMPING FOR THE MISSING SNPS
grep -w -A10 'everymissingsnip' rsqSNPseduyears

#11.EXTRACT THE PROXY SNPS THAT WERE FOUND FOR ANY OF THE MISSING SNPs
grep -A45 'rs1396967' rsqSNPseduyears | grep -v INDEX | cut -d ' ' -f2,4,7 > temp
awk '$2 > 0.96' temp | cut -d ' ' -f1 > rs1396967rsqSNPseduyears
grep -A19 'rs111321694' rsqSNPseduyears | grep -v INDEX | cut -d ' ' -f2,4,7 > temp
awk '$2 > 0.975 && $3 < .0000002' temp | cut -d ' ' -f1 > rs111321694rsqSNPseduyears
grep -A2026 'rs192818565' rsqSNPseduyears | grep -v INDEX | cut -d ' ' -f2,4,7 > temp
awk '$2 > 0.93 && $3 < .00000005' temp | cut -d ' ' -f1 > rs192818565rsqSNPseduyears
grep -A23 'rs28792186' rsqSNPseduyears | grep -v INDEX | cut -d ' ' -f2,4,7 > temp
awk '$2 > 0.97' temp | cut -d ' ' -f1 > rs28792186rsqSNPseduyears
grep -A4 'rs6774721' rsqSNPseduyears | grep -v INDEX | cut -d ' ' -f2 > rs6774721rsqSNPseduyears
grep -A6 'rs7776010' rsqSNPseduyears | grep -v INDEX | cut -d ' ' -f2 > rs7776010rsqSNPseduyears
grep -A4 'rs12534506' rsqSNPseduyears | grep -v INDEX | cut -d ' ' -f2 > rs12534506rsqSNPseduyears
grep -A1 'rs1106761' rsqSNPseduyears | grep -v INDEX | cut -d ' ' -f2 > rs1106761rsqSNPseduyears

#12.EXTRAXT THIS NEW PROXY SNPs FROM THE PD GWAS DATASET, IF THERE ARE FOUND
grep -w -f rs1396967rsqSNPseduyears /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs1396967rsqSNPseduyearsPD
grep -w -f rs111321694rsqSNPseduyears /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs111321694rsqSNPseduyearsPD
grep -w -f rs192818565rsqSNPseduyears /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs192818565rsqSNPseduyears
grep -w -f rs28792186rsqSNPseduyears /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs28792186rsqSNPseduyears
grep -w -f rs6774721rsqSNPseduyears /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs6774721rsqSNPseduyears
grep -w -f rs7776010rsqSNPseduyears /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs7776010rsqSNPseduyears
grep -w -f rs12534506rsqSNPseduyears /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs12534506rsqSNPseduyears
grep -w -f rs1106761rsqSNPseduyears /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs1106761rsqSNPseduyears

#13.JOIN ALL INDEX SNPS FOUND ON THE PD GWAS DATASET
cat eduyearsSNPinPDdataV rs1396967rsqSNPseduyearsPD rs111321694rsqSNPseduyearsPD rs192818565rsqSNPseduyears rs28792186rsqSNPseduyears rs6774721rsqSNPseduyears rs7776010rsqSNPseduyears rs12534506rsqSNPseduyears rs1106761rsqSNPseduyears > fulleduyearsSNPinPDdataV

#14.EXTRACT THE INFO FOR CALCULATING Rsquare AND F-STATISTIC FROM INDEX SNPS ALSO PRESENT IN THE PD GWAS DATASET 
echo "MarkerName CHR POS A1 A2 EAF Beta SE Pval" | sed -e 's/ /\t/g'  > header
awk '{ if (NR!=1) print $1}' fulleduyearsSNPinPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/education/EduYears_Main.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/eduyearsforR2calcV
