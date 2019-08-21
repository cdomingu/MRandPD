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

#7. AS PD NALLS DATASET DOESN'T HAVE AN rsID, EXTRACT THE CHR:BP POSITION
cd /data/neurogen/MRandPD/Results/R2andFstatisticprePD/
awk '{ if (NR!=1) print $1,$2}' drinksperweekforR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > drinkspweekchrpos

echo "SNP rsID" > header2
awk '{ if (NR!=1) print $3}' drinksperweekforR2calc > drinkspweekrsID
paste drinkspweekchrpos drinkspweekrsID | sort > temp

#8. EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET
echo "SNP A1 A2 freq b se p N_cases N_controls" > header
grep -w -f drinkspweekchrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/drinkspweekSNPinNallsPDdataV

#9. AD THE rsID
join drinkspweekSNPinNallsPDdataV /data/neurogen/MRandPD/Results/R2andFstatisticprePD/drinkspweekchrposrsID | sed -e 's/ /\t/g' > temp
mv temp drinkspweekSNPinNallsPDdataV

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

#7. AS PD NALLS DATASET DOESN'T HAVE AN rsID, EXTRACT THE CHR:BP POSITION
cd /data/neurogen/MRandPD/Results/R2andFstatisticprePD/
awk '{ if (NR!=1) print $2,$3}' eduyearsforR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > eduyearschrpos

echo "SNP rsID" > header2
awk '{ if (NR!=1) print $3}' drinksperweekforR2calc > drinkspweekrsID
paste drinkspweekchrpos drinkspweekrsID | sort > temp

#8. EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET
echo "SNP A1 A2 freq b se p N_cases N_controls" > header
grep -w -f eduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/eduyearsSNPinNallsPDdataV

#9. GET MISSING VARIANTS (INDEX SNPS THAT WEREN'T FOUND ON PD GWAS DATASET)
cd /data/neurogen/MRandPD/Results/HarmonizatioNalls2019/
awk '{ if (NR!=1) print $1}' eduyearsSNPinNallsPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/R2andFstatisticprePD/eduyearschrpos | sort  > missingeduyearsSNPinNallsPD

grep -f missingeduyearsSNPinNallsPD /data/neurogen/MRandPD/Results/R2andFstatisticprePD/eduyearsposrsID > temp
mv temp missingeduyearsSNPinNallsPD

#10.IDENTIFY PROXY SNPS FROM MISSING VARIANT (THE ONES WITH Rsquare > 0.9 WITH INDEX SNP)
cd /data/neurogen/MRandPD/Results/Clumping/EduyearsClumppvalmin8Verbose/
cat clumped.eduyears.chr*.clumped > temp
sed -e 's/     / /g' temp| sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > prueba
awk -F'[ ]' '$4>0.9 || $5==1' prueba | grep 'rs' | grep -v "^\s[0-9]" > rsqSNPseduyears
awk -F'[ ]' '$4>0.7 || $5==1' prueba | grep 'rs' | grep -v "^\s[0-9]" > rsq7SNPseduyears
awk -F'[ ]' '$4>0.6 || $5==1' prueba | grep 'rs' | grep -v "^\s[0-9]" > rsq5SNPseduyears

#11."MANUALLY" DETERMINE HOW MANY PROXY SNPS WERE FOUND DURING THE CLUMPING FOR THE MISSING SNPS
grep -w -A10 'everymissingsnip' rsqSNPseduyears

grep -w -A20 'rs11222416' rsq7SNPseduyears | grep -v INDEX | cut -d ' ' -f2,4,7 > temp
awk '$2 > 0.75 && $3 < .0000004' temp | cut -d ' ' -f1 > rs11222416rsqSNPseduyears
grep -w -A37 'rs7948975' rsq7SNPseduyears | grep -v INDEX | cut -d ' ' -f2,4,7 > temp
awk '$3 < .00000008' temp | cut -d ' ' -f1 > rs7948975rsqSNPseduyears
grep -w -A120 'rs12900061' rsqSNPseduyears | grep -v INDEX | cut -d ' ' -f2,4,7 > temp
awk '$2 > 0.99 && $3 < .000000006' temp | cut -d ' ' -f1 > rs12900061rsqSNPseduyears
grep -w -A24 'rs28420834' rsq7SNPseduyears | grep -v INDEX | cut -d ' ' -f2,4,7 > temp
awk '$2 > 0.72 && $3 < .0000002' temp | cut -d ' ' -f1 > rs28420834rsqSNPseduyears
grep -w -A13 'rs11726992' rsq7SNPseduyears | grep -v INDEX | cut -d ' ' -f2,4,7 > temp
awk '$2 > 0.74 && $3 < .000006' temp | cut -d ' ' -f1 > rs11726992rsqSNPseduyears
grep -w -A5 'rs9792504' rsqSNPseduyears | grep -v INDEX | cut -d ' ' -f2 > rs9792504rsqSNPseduyears
grep -w -A8 'rs7033137' rsq6SNPseduyears | grep -v INDEX | cut -d ' ' -f2,4,7 > temp
awk '$3 < .0000006' temp | cut -d ' ' -f1 > rs7033137rsqSNPseduyears

#12.EXTRACT THE PROXY SNPS THAT WERE FOUND FOR ANY OF THE MISSING SNPs
cut -f1,2,3 /data/neurogen/MRandPD/GWASsummaryStats/education/EduYears_Main.txt > EduYearschrposID
grep -w -f rs11222416rsqSNPseduyears EduYearschrposID > rs11222416proxySNPeduyearschrposID
cut -f2,3 rs11222416proxySNPeduyearschrposID | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > rs11222416SNPeduyearschrpos
grep -w -f rs7948975rsqSNPseduyears EduYearschrposID > rs7948975proxySNPeduyearschrposID
cut -f2,3 rs7948975proxySNPeduyearschrposID | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > rs7948975SNPeduyearschrpos
grep -w -f rs12900061rsqSNPseduyears EduYearschrposID > rs12900061proxySNPeduyearschrposID
cut -f2,3 rs12900061proxySNPeduyearschrposID | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > rs12900061SNPeduyearschrpos
grep -w -f rs28420834rsqSNPseduyears EduYearschrposID > rs28420834proxySNPeduyearschrposID
cut -f2,3 rs28420834proxySNPeduyearschrposID | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > rs28420834SNPeduyearschrpos
grep -w -f rs11726992rsqSNPseduyears EduYearschrposID > rs11726992proxySNPeduyearschrposID
cut -f2,3 rs11726992proxySNPeduyearschrposID | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > rs11726992SNPeduyearschrpos
grep -w -f rs9792504rsqSNPseduyears EduYearschrposID > rs9792504proxySNPeduyearschrposID
cut -f2,3 rs9792504proxySNPeduyearschrposID | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > rs9792504SNPeduyearschrpos
grep -w -f rs7033137rsqSNPseduyears EduYearschrposID > rs7033137proxySNPeduyearschrposID
cut -f2,3 rs7033137proxySNPeduyearschrposID | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > rs7033137SNPeduyearschrpos
grep -w -f rs1396967rsqSNPseduyears EduYearschrposID > rs1396967proxySNPeduyearschrposID
cut -f2,3 rs1396967proxySNPeduyearschrposID | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > rs1396967SNPeduyearschrpos
grep -w -f rs111321694rsqSNPseduyears EduYearschrposID > rs111321694proxySNPeduyearschrposID
cut -f2,3 rs111321694proxySNPeduyearschrposID | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > rs111321694SNPeduyearschrpos
grep -w -f rs192818565rsqSNPseduyears EduYearschrposID > rs192818565proxySNPeduyearschrposID
cut -f2,3 rs192818565proxySNPeduyearschrposID | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > rs192818565SNPeduyearschrpos
grep -w -f rs71537331rsqSNPseduyears EduYearschrposID > rs71537331proxySNPeduyearschrposID
cut -f2,3 rs71537331proxySNPeduyearschrposID | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > rs71537331SNPeduyearschrpos
grep -w -f rs7776010rsqSNPseduyears EduYearschrposID > rs7776010proxySNPeduyearschrposID
cut -f2,3 rs7776010proxySNPeduyearschrposID | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > rs7776010SNPeduyearschrpos
grep -w -f rs28792186rsqSNPseduyears EduYearschrposID > rs28792186proxySNPeduyearschrposID
cut -f2,3 rs28792186proxySNPeduyearschrposID | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > rs28792186SNPeduyearschrpos

#13.EXTRAXT THIS NEW PROXY SNPs FROM THE PD GWAS DATASET, IF THERE ARE FOUND
grep -w -f rs11222416SNPeduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyears/rs11222416proxSNPseduyearsPD
grep -w -f rs7948975SNPeduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyears/rs7948975proxSNPseduyearsPD
grep -w -f rs12900061SNPeduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyears/rs12900061proxSNPseduyearsPD
grep -w -f rs28420834SNPeduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyears/rs28420834proxSNPseduyearsPD
grep -w -f rs11726992SNPeduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyears/rs11726992proxSNPseduyearsPD
grep -w -f rs9792504SNPeduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyears/rs9792504proxSNPseduyearsPD
grep -w -f rs7033137SNPeduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyears/rs7033137proxSNPseduyearsPD
grep -w -f rs1396967SNPeduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyears/rs1396967proxSNPseduyearsPDPD
grep -w -f rs111321694SNPeduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyears/rs111321694proxSNPseduyearsPD
grep -w -f rs192818565SNPeduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyears/rs192818565proxSNPseduyearsPD
grep -w -f rs71537331SNPeduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyears/rs71537331proxSNPseduyearsPD
grep -w -f rs7776010SNPeduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyears/rs7776010proxSNPseduyearsPD
grep -w -f rs28792186SNPeduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyearsrs28792186proxSNPseduyearsPD

#14.JOIN ALL INDEX SNPS FOUND ON THE PD GWAS DATASET
cat eduyearsSNPinNallsPDdataV /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ProxySNPseduyears/*proxSNPseduyears* | sort | uniq > fulleduyearsSNPinNallsPDdataV 

#15. ADD THE rsID TO THE INFO OF THE INDEX SNPS ON PD DATA
cd /data/neurogen/MRandPD/GWASsummaryStats/education/
echo "SNP rsID" > header2
awk '{ if (NR!=1) print $2,$3}' EduYears_Main.txt > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > GWASeduyearschrpos
awk '{ if (NR!=1) print $1}' EduYears_Main.txt > GWASeduyearsID
paste GWASeduyearschrpos GWASeduyearsID | sort > temp
cat header2 temp | sed -e 's/ /\t/g' > GWASeduyearschrposrsID

cd /data/neurogen/MRandPD/Results/HarmonizatioNalls2019/
awk '{ if (NR!=1) print $1}' fulleduyearsSNPinNallsPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/education/GWASeduyearschrposrsID | sort > fulleduyearsSNPinNallsPDdataVchrposID
join fulleduyearsSNPinNallsPDdataV fulleduyearsSNPinNallsPDdataVchrposID > temp
cat header3 temp | sed -e 's/ /\t/g' > fulleduyearsSNPinNallsPDdataV

#16.EXTRACT THE INFO FOR CALCULATING Rsquare AND F-STATISTIC FROM INDEX SNPS ALSO PRESENT IN THE PD GWAS DATASET 
echo "MarkerName CHR POS A1 A2 EAF Beta SE Pval" | sed -e 's/ /\t/g'  > header
awk '{ if (NR!=1) print $10}' fulleduyearsSNPinNallsPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/education/EduYears_Main.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/eduyearsforR2calcVNalls
