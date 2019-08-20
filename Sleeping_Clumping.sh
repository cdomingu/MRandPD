##This code is for obtainig the index SNPs (IV) by clumping from SLEEPING GWAS data
##We need to clump SNPs to ensure that all the index SNPs are independent frome one another
##After that, we will extract those index SNPs data from the PD GWAS, an calculate the Rsquare an F-statistic from the instruments
##We are using the summary statistics from Lane,2017 https://doi.org/10.1038/ng.3749

# 2. PREPARE GWAS FILES FOR SNP CLUMPING (in PLINK)
# EXTRACTION OF COLUMNS (SNP EA NEA P BETA) from the GWAS results
#Issue with chromosome 23? in sleeping disturbances data
#get rid of the extra info following the snp rs number, select just snp bp pval columns and re run test
cd /data/neurogen/MRandPD/GWASsummaryStats/sleep/
echo "SNP BP P" | sed -e 's/ /\t/g'  > header
awk '{ if (NR!=1) print $1,$3,$9}' Application6818_SleepTrait_data_Lane_RESULTS_DATA_EXCESSIVE_DAYTIME_SLEEPINESS_all.txt | sed -e 's/ /\t/g' > temp
grep ':' temp  > issue
grep -F -x -v -f issue temp > data
sed -e 's/:/\t/g' issue | cut -f1,5,6 > temp
cat header data temp > /data/neurogen/MRandPD/Results/Clumping/Data4Clump/excessdaysleepallwo23
awk '{ if (NR!=1) print $1,$3,$9}' Application6818_SleepTrait_data_Lane_RESULTS_DATA_INSOMNIA_SYMPTOMS_all.txt | sed -e 's/ /\t/g' > temp
grep ':' temp  > issue
grep -F -x -v -f issue temp > data
sed -e 's/:/\t/g' issue | cut -f1,5,6 > temp
cat header data temp > /data/neurogen/MRandPD/Results/Clumping/Data4Clump/insomniallwo23
awk '{ if (NR!=1) print $1,$3,$9}' Application6818_SleepTrait_data_Lane_RESULTS_DATA_SLEEP_DURATION_all.txt | sed -e 's/ /\t/g' > temp
grep ':' temp  > issue
grep -F -x -v -f issue temp > data
sed -e 's/:/\t/g' issue | cut -f1,5,6 > temp
cat header data temp > /data/neurogen/MRandPD/Results/Clumping/Data4Clump/sleepdurallwo23

# 3.1 COLLATE ALL SNPS FROM GWAS
# For sleeping GWAS
cd /data/neurogen/MRandPD/Results/Clumping/Data4Clump/
for i in excessdaysleepallwo23 insomniallwo23 sleepdurallwo23 
do echo $i
awk '{ if (NR != 1) print $1}' $i >> allSnpRsGWASsleepwo23
done
# 3.3 SORT RS NUMBERS AND REMOVE DUPLICATED SNPS
sort allSnpRsGWASsleepwo23 | uniq > allSnpRsUniquesleepwo23
# 3.4 GET RID OF INDELS (THEY DO NOT HAVE AN RS NUMBER)
grep 'rs' allSnpRsUniquesleepwo23 > allSnpRsUniqueNoIndelssleepwo23

# We will be ussing the 1000Genomes data in Plink format downloaded from: http://www.cog-genomics.org/plink/1.9/resources for clumping

# 4 CREATE .lsf FILES TO RUN IN CLUSTER QUEUE
#Parameters for clumping: --clump-p1 .0000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose
#for identifying both genomewide significant and suggestive index snps which are independent form one another
#this mean that intex snps will reach significance (.0000005) and either have a r2=0.001 with the next index SNP or be 10000 kb appart
#--clump-verbose to ger the Rsquare values of the clumped SNPs with the index SNP
for i in excessdaysleepallwo23 insomniallwo23 sleepdurallwo23
do
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
cp /data/neurogen/MRandPD/Scripts/test.lsf /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/SleepClumpScrptpvalmin7Verbose/ 
mv /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/SleepClumpScrptpvalmin7Verbose/test.lsf /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/SleepClumpScrptpvalmin7Verbose/"$i".clump_chr"$j".lsf
echo "module load plink/1.90b3" >> /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/SleepClumpScrptpvalmin7Verbose/"$i".clump_chr"$j".lsf
echo "plink --bfile /data/neurogen/MRandPD/1000GenomsPlinkFormat/1kg_phase1_chr"$j" --extract /data/neurogen/MRandPD/Results/Clumping/Data4Clump/allSnpRsUniqueNoIndelssleepwo23 --clump /data/neurogen/MRandPD/Results/Clumping/Data4Clump/"$i" --clump-p1 .0000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose --out /data/neurogen/MRandPD/Results/Clumping/SleepClumpvalmin7Verbose/clumped."$i".chr"$j" " >> /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/SleepClumpScrptpvalmin7Verbose/"$i".clump_chr"$j".lsf
done
done

# SUBMIT THEM TO THE QUEUE
for i in excessdaysleepallwo23 insomniallwo23 sleepdurallwo23
do
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 
do
bsub -q normal < "$i".clump_chr"$j".lsf
done
done

#5.EXTRACT THE LIST OF INDEXSNPS FROM CLUMPING RESULTS
cd /data/neurogen/MRandPD/Results/Clumping/SleepClumpvalmin7Verbose/
cat clumped.excessdaysleepallwo23.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPSexcessdaysleepV
cat clumped.insomniallwo23.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPSinsomniaV
cat clumped.sleepdurallwo23.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPSsleepdurV

#6.EXTRACT EXPOSURE GWAS INFO (specially BETA, MAF,SE) FOR CALCUATING PHENOTYPIC VARIANCE EXPALINED (Rsquare) AND F-STATISTIC
echo "SNP CHR BP A1 A2 MAF BETA SE P N" | sed -e 's/ /\t/g'  > header
grep -w -f indexSNPSexcessdaysleepV /data/neurogen/MRandPD/GWASsummaryStats/sleep/Application6818_SleepTrait_data_Lane_RESULTS_DATA_EXCESSIVE_DAYTIME_SLEEPINESS_all.txt | sort -n -k 2 | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticprePD/excessdaysleepforR2calc
grep -w -f indexSNPSinsomniaV /data/neurogen/MRandPD/GWASsummaryStats/sleep/Application6818_SleepTrait_data_Lane_RESULTS_DATA_INSOMNIA_SYMPTOMS_all.txt | sort -n -k 2 | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticprePD/insomniaforR2calc
grep -w -f indexSNPSsleepdurV /data/neurogen/MRandPD/GWASsummaryStats/sleep/Application6818_SleepTrait_data_Lane_RESULTS_DATA_SLEEP_DURATION_all.txt | sort -n -k 2 | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticprePD/sleepdurforR2calc

#7. AS PD DATASET DOS NOT HAVE rsID BUT CHR:BP VALUES, EXTRACT THE CHR:POS OF OUR INDEX SNPS
awk '{ if (NR!=1) print $2,$3}' excessdaysleepforR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > excessdaysleepchrpos
awk '{ if (NR!=1) print $2,$3}' sleepdurforR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > sleepdurchrpos
awk '{ if (NR!=1) print $2,$3}' insomniaforR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > insomniachrpos
awk '{ if (NR!=1) print $2,$3}' eduyearsforR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > eduyearschrpos

#8.EXTRACT OUR INDEX SNPS CHR:BP FROM PD DATASET
echo "SNP A1 A2 freq b se p N_cases N_controls" > header
grep -w -f excessdaysleepchrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab| sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/excessdaysleepSNPinNallsPDdataV
grep -w -f insomniachrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/insomniaSNPinNallsPDdataV
grep -w -f sleepdurchrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/sleepdurSNPinNallsPDdataV

#9. MAKE A rsID-CH:POS CONVERTION TABLE OF OUR INDEX SNPs
awk '{ if (NR!=1) print $1}' excessdaysleepforR2calc > excessdaysleeprsID
paste excessdaysleepchrpos excessdaysleeprsID | sort > temp
cat header2 temp | sed -e 's/ /\t/g' > excessdaysleepposrsID
awk '{ if (NR!=1) print $1}' sleepdurforR2calc > sleepdurrsID
paste sleepdurchrpos sleepdurrsID | sort > temp
cat header2 temp | sed -e 's/ /\t/g' > sleepdurposrsID
awk '{ if (NR!=1) print $1}' insomniaforR2calc > insomniarsID
paste insomniachrpos insomniarsID | sort > temp
cat header2 temp | sed -e 's/ /\t/g' > insomniaposrsID

#10.GET MISSING VARIANTS (INDEX SNPS THAT WEREN'T FOUND ON PD GWAS DATASET)
cd /data/neurogen/MRandPD/Results/HarmonizationNalls2019
awk '{ if (NR!=1) print $1}' excessdaysleepSNPinNallsPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/R2andFstatisticprePD/excessdaysleepchrpos | sort > missingexcessdaysleepSNPinNallsPD
awk '{ if (NR!=1) print $1}' insomniaSNPinNallsPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/R2andFstatisticprePD/insomniachrpos | sort > missinginsomniaSNPinNallsPD
awk '{ if (NR!=1) print $1}' sleepdurSNPinNallsPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/R2andFstatisticprePD/sleepdurchrpos | sort > missingsleepdurSNPinNallsPD

#11.IDENTIFY PROXY SNPS FROM MISSING VARIANT (THE ONES WITH Rsquare > 0.6 WITH INDEX SNP)
cd /data/neurogen/MRandPD/Results/Clumping/SleepClumpvalmin7Verbose/
cat clumped.excessdaysleepallwo23.chr*.clumped > temp
sed -e 's/     / /g' temp| sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > test
awk -F'[ ]' '$4>0.9 || $5==1' prueba | grep 'rs' | grep -v "^\s[0-9]" > rsqSNPsexcessdaysleep
awk -F'[ ]' '$4>0.8 || $5==1' prueba | grep 'rs' | grep -v "^\s[0-9]" > rsq8SNPsexcessdaysleep
awk -F'[ ]' '$4>0.7 || $5==1' prueba | grep 'rs' | grep -v "^\s[0-9]" > rsq7SNPsexcessdaysleep
awk -F'[ ]' '$4>0.6 || $5==1' prueba | grep 'rs' | grep -v "^\s[0-9]" > rsq6SNPsexcessdaysleep

cat clumped.insomniallwo23.chr*.clumped > test
sed -e 's/     / /g' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > test
awk -F'[ ]' '$4>0.9 || $5==1' test | grep 'rs' | grep -v "^\s[0-9]" > rsqSNPsinsomnia
awk -F'[ ]' '$4>0.6 || $5==1' prueba | grep 'rs' | grep -v "^\s[0-9]" > rsq6SNPsinsomnia
#after manually checking, no new snps for snp duration

cat clumped.sleepdurallwo23.chr*.clumped > temp
sed -e 's/     / /g' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > test
awk -F'[ ]' '$4>0.6 || $5==1' test | grep 'rs' | grep -v "^\s[0-9]" > rsqSNPsleepdurall
#after manually checking, no new snps for snp duration

#12."MANUALLY" DETERMINE HOW MANY PROXY SNPS WERE FOUND DURING THE CLUMPING FOR THE MISSING SNPS
grep -w -A10 'everymissingsnp' rsqSNPsexcessdaysleep
#no proxy snps were found for the next variants:
#rs147917224
#rs148846329
#rs188800342
#rs184375653
#rs2681335
grep -w A10 'everymissingsnp' rsqSNPsinsomnia

#13.EXTRACT THE PROXY SNPS THAT WERE FOUND FOR ANY OF THE MISSING SNPs
grep -A9 'rs11851923' rsq6SNPsexcessdaysleep | grep -v INDEX | cut -d ' ' -f2 > rs11851923RsqSNPexcessdaysleep
grep -A1 'rs116450109' rsq6SNPsexcessdaysleep | grep -v INDEX | cut -d ' ' -f2 > rs116450109RsqSNPexcessdaysleep
grep -A5 'rs115320831' rsq7SNPsexcessdaysleep | grep -v INDEX | cut -d ' ' -f2 > rs115320831RsqSNPexcessdaysleep
grep -A2 'rs35309287' rsq7SNPsexcessdaysleep | grep -v INDEX | cut -d ' ' -f2 > rs35309287RsqSNPexcessdaysleep
grep -A1 'rs143088939' rsq8SNPsexcessdaysleep | grep -v INDEX | cut -d ' ' -f2 > rs143088939RsqSNPexcessdaysleep

cut -f1,2,3  -d ' ' /data/neurogen/MRandPD/GWASsummaryStats/sleep/Application6818_SleepTrait_data_Lane_RESULTS_DATA_EXCESSIVE_DAYTIME_SLEEPINESS_all.txt > ESCESSDAYSLEEPchrposID
grep -w -f rs115320831RsqSNPexcessdaysleep ESCESSDAYSLEEPchrposID > rs115320831proxySNPexcessdaysleepID
cut -f2,3 -d ' ' rs115320831proxySNPexcessdaysleepID | sed -e 's/^/chr/g' | sed -e 's/ /:/g' > rs115320831SNPexcessdaysleepchrpos
#grep -w -f rs11851923RsqSNPexcessdaysleep ESCESSDAYSLEEPchrposID > rs11851923proxySNPexcessdaysleepID
#cut -f2,3 -d ' ' rs11851923proxySNPexcessdaysleepID | sed -e 's/^/chr/g' | sed -e 's/ /:/g' > rs11851923SNPexcessdaysleepchrpos
#grep -w -f rs116450109RsqSNPexcessdaysleep ESCESSDAYSLEEPchrposID > rs116450109proxySNPexcessdaysleepID
#cut -f2,3 -d ' ' rs116450109proxySNPexcessdaysleepID | sed -e 's/^/chr/g' | sed -e 's/ /:/g' > rs116450109SNPexcessdaysleepchrpos
grep -w -f rs35309287RsqSNPexcessdaysleep ESCESSDAYSLEEPchrposID > rs35309287proxySNPexcessdaysleepID
cut -f2,3 -d ' ' rs35309287proxySNPexcessdaysleepID | sed -e 's/^/chr/g' | sed -e 's/ /:/g' > rs35309287SNPexcessdaysleepchrpos
#grep -w -f rs143088939RsqSNPexcessdaysleep ESCESSDAYSLEEPchrposID > rs143088939proxySNPexcessdaysleepID
#cut -f2,3 -d ' ' rs143088939proxySNPexcessdaysleepID | sed -e 's/^/chr/g' | sed -e 's/ /:/g' > rs143088939SNPexcessdaysleepchrpos

#14. LOOK FOR THEM ON THE PD DATASET
#grep -w -f rs11851923SNPexcessdaysleepchrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/rs11851923proxSNPexcessdaysleepPD
#grep -w -f rs116450109SNPexcessdaysleepchrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/rs116450109proxSNPexcessdaysleepPD
grep -w -f rs35309287SNPexcessdaysleepchrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/rs35309287proxSNPexcessdaysleepPD
#grep -w -f rs143088939SNPexcessdaysleepchrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/rs143088939proxSNPexcessdaysleepPD
grep -w -f rs115320831SNPexcessdaysleepchrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/rs11532083proxSNPexcessdaysleepPD

#15. JOIN THE RESULTS OF THE PROXY SNPS WITH THE REST OF THE INDEX SNPS, AND THEN AD THE rsID
cd /data/neurogen/MRandPD/Results/HarmonizationNalls2019
cat excessdaysleepSNPinNallsPDdataV rs35309287proxSNPexcessdaysleepPD rs11532083proxSNPexcessdaysleepPD | sort | uniq > fullexcessdaysleepSNPinNallsPDdataV

join insomniaSNPinNallsPDdataV /data/neurogen/MRandPD/Results/R2andFstatisticprePD/insomniaposrsID | sed -e 's/ /\t/g' > temp
mv temp insomniaSNPinNallsPDdataV
join sleepdurSNPinNallsPDdataV /data/neurogen/MRandPD/Results/R2andFstatisticprePD/sleepdurposrsID | sed -e 's/ /\t/g' > temp
mv temp sleepdurSNPinNallsPDdataV

echo "SNP A1 A2 freq b se p N_cases N_controls rsID" | sed -e 's/ /\t/g' > header3
awk '{ if (NR!=1) print $1}' fullexcessdaysleepSNPinNallsPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/DBsnpchrposID | sort > fullexcessdaysleepSNPinNallsPDdataVchrposID
join fullexcessdaysleepSNPinNallsPDdataV fullexcessdaysleepSNPinNallsPDdataVchrposID > temp
cat header3 temp | sed -e 's/ /\t/g' > fullexcessdaysleepSNPinNallsPDdataV

#16. THEN WE WILL EXTRACT THE INFORMATION FOR THIS NEW SET TO SNPs TO CALCULATE THEIR Rsq AND F-STATISTIC
echo "SNP CHR BP A1 A2 MAF BETA SE P N" | sed -e 's/ /\t/g'  > header
awk '{ if (NR!=1) print $10}' fullexcessdaysleepSNPinNallsPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/sleep/Application6818_SleepTrait_data_Lane_RESULTS_DATA_EXCESSIVE_DAYTIME_SLEEPINESS_all.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/fullexcessdaysleepforR2calcVNalls
awk '{ if (NR!=1) print $10}' insomniaSNPinNallsPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/sleep/Application6818_SleepTrait_data_Lane_RESULTS_DATA_INSOMNIA_SYMPTOMS_all.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/insomniaforR2calcVNalls
awk '{ if (NR!=1) print $10}' sleepdurSNPinNallsPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/education/EduYears_Main.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/sleepdurforR2calcVNalls

