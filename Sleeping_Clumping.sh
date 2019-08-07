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

#7. EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET
echo "MarkerName Allele1 Allele2 Effect StdErr P-value NStudies HetISq" | sed -e 's/ /\t/g'  > header
grep -w -f indexSNPSexcessdaysleepV /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/excessdaysleepSNPinPDdataV
grep -w -f indexSNPSinsomniaV /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/insomniaSNPinPDdataV
grep -w -f indexSNPSsleepdurV /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/sleepdurSNPinPDdataV

#8.GET MISSING VARIANTS (INDEX SNPS THAT WEREN'T FOUND ON PD GWAS DATASET)
cd /data/neurogen/MRandPD/Results/HarmonizationVerbose/
awk '{ if (NR!=1) print $1}' excessdaysleepSNPinPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/Clumping/SleepClumpvalmin7Verbose/indexSNPSexcessdaysleepV > missingindexSNPexcessdaysleepV
awk '{ if (NR!=1) print $1}' insomniaSNPinPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/Clumping/SleepClumpvalmin7Verbose/indexSNPSinsomniaV > missingindexSNPinsomniaV
awk '{ if (NR!=1) print $1}' sleepdurSNPinPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/Clumping/SleepClumpvalmin7Verbose/indexSNPSsleepdurV > missingindexSNPsleepdur

#9.IDENTIFY PROXY SNPS FROM MISSING VARIANT (THE ONES WITH Rsquare > 0.9 WITH INDEX SNP)
cd /data/neurogen/MRandPD/Results/Clumping/SleepClumpvalmin7Verbose/
cat clumped.excessdaysleepallwo23.chr*.clumped > temp
sed -e 's/     / /g' temp| sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > test
awk -F'[ ]' '$4>0.9 || $5==1' prueba | grep 'rs' | grep -v "^\s[0-9]" > rsqSNPsexcessdaysleep
cat clumped.insomniallwo23.chr*.clumped > test
sed -e 's/     / /g' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > test
awk -F'[ ]' '$4>0.9 || $5==1' test | grep 'rs' | grep -v "^\s[0-9]" > rsqSNPsinsomnia
#after manually checking, no new snps
cat clumped.sleepdurallwo23.chr*.clumped > temp
sed -e 's/     / /g' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > test
awk -F'[ ]' '$4>0.9 || $5==1' test | grep 'rs' | grep -v "^\s[0-9]" > rsqSNPsleepdurall
#after manually checking, no new snps

#10."MANUALLY" DETERMINE HOW MANY PROXY SNPS WERE FOUND DURING THE CLUMPING FOR THE MISSING SNPS
grep -w -A6 'everymissingsnp' rsqSNPsexcessdaysleep

#11.EXTRACT THE PROXY SNPS THAT WERE FOUND FOR ANY OF THE MISSING SNPs
grep -A2 'rs180885632' rsqSNPsexcessdaysleep | grep -v INDEX | cut -d ' ' -f2 > rs180885632RsqSNPexcessdaysleep
grep -A1 'rs11851923' rsqSNPsexcessdaysleep | grep -v INDEX | cut -d ' ' -f2 > rs11851923RsqSNPexcessdaysleep
grep -A2 'rs192315283' rsqSNPsexcessdaysleep | grep -v INDEX | cut -d ' ' -f2 > rs192315283RsqSNPexcessdaysleep
grep -A2 'rs138066714' rsqSNPsexcessdaysleep | grep -v INDEX | cut -d ' ' -f2 > rs138066714RsqSNPexcessdaysleep
grep -A1 'rs35309287' rsqSNPsexcessdaysleep | grep -v INDEX | cut -d ' ' -f2 > rs35309287RsqSNPexcessdaysleep
grep -A1 'rs71564440' rsqSNPsexcessdaysleep | grep -v INDEX | cut -d ' ' -f2 > rs71564440RsqSNPexcessdaysleep
grep -A1 'rs58357132' rsqSNPsexcessdaysleep | grep -v INDEX | cut -d ' ' -f2 > rs58357132RsqSNPexcessdaysleep

#12.EXTRAXT THIS NEW PROXY SNPs FROM THE PD GWAS DATASET, IF THERE ARE FOUND
#grep -w -f rs180885632RsqSNPexcessdaysleep /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs180885632RsqSNPexcessdaysleepPD
#grep -w -f rs11851923RsqSNPexcessdaysleep /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs11851923RsqSNPexcessdaysleepPD
#grep -w -f rs192315283RsqSNPexcessdaysleep /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs192315283RsqSNPexcessdaysleepPD
#grep -w -f rs138066714RsqSNPexcessdaysleep /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs138066714RsqSNPexcessdaysleepPD
#grep -w -f rs35309287RsqSNPexcessdaysleep /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs35309287RsqSNPexcessdaysleepPD
grep -w -f rs71564440RsqSNPexcessdaysleep /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs71564440RsqSNPexcessdaysleepPD
#grep -w -f rs58357132RsqSNPexcessdaysleep /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs58357132RsqSNPexcessdaysleepPD

#13.JOIN ALL INDEX SNPS FOUND ON THE PD GWAS DATASET
cat excessdaysleepSNPinPDdataV rs71564440RsqSNPexcessdaysleepPD > fullexcessdaysleepSNPinPDV

#14.EXTRACT THE INFO FOR CALCULATING Rsquare AND F-STATISTIC FROM INDEX SNPS ALSO PRESENT IN THE PD GWAS DATASET 
echo "SNP CHR BP A1 A2 MAF BETA SE P N" | sed -e 's/ /\t/g'  > header
awk '{ if (NR!=1) print $1}' fullexcessdaysleepSNPinPDV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/sleep/Application6818_SleepTrait_data_Lane_RESULTS_DATA_EXCESSIVE_DAYTIME_SLEEPINESS_all.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/excessdaysleepforR2calcV
awk '{ if (NR!=1) print $1}' insomniaSNPinPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/sleep/Application6818_SleepTrait_data_Lane_RESULTS_DATA_INSOMNIA_SYMPTOMS_all.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/insomniaforR2calcV
awk '{ if (NR!=1) print $1}' sleepdurSNPinPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/sleep/Application6818_SleepTrait_data_Lane_RESULTS_DATA_SLEEP_DURATION_all.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/sleepdurforR2calcV
