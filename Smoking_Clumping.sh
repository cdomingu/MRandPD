##This code is for obtainig the index SNPs (IV) by clumping from smoking GWAS data
##We need to clump SNPs to ensure that all the index SNPs are independent frome one another
##Ater that, we will extract those index SNPs data from the PD GWAS, an calculate the Rsquare an F-statistic from the instruments
##We are using the summary statistics from Liu,2019 downloaded from https://conservancy.umn.edu/handle/11299/201564

#2. PREPARE GWAS FILES FOR SNP CLUMPING (in PLINK)
# EXTRACTION OF COLUMNS (SNP CHR BP P) from the GWAS results
cd /data/neurogen/MRandPD/GWASsummaryStats/Smoking2019/
echo "SNP CHR BP P" | sed -e 's/ /\t/g'  > header
awk '{ if (NR!=1) print $3,$1,$2,$8}' CigarettesPerDay.txt | sed -e 's/ /\t/g' > temp
cat header temp > /data/neurogen/MRandPD/Results/Clumping/Data4Clump/cigarettespday
awk '{ if (NR!=1) print $3,$1,$2,$8}' EverSmoking.txt | sed -e 's/ /\t/g' > temp
cat header temp > /data/neurogen/MRandPD/Results/Clumping/Data4Clump/eversmoking
awk '{ if (NR!=1) print $3,$1,$2,$8}' SmokingStatus.txt | sed -e 's/ /\t/g' > temp
cat header temp > /data/neurogen/MRandPD/Results/Clumping/Data4Clump/smokingstatus
awk '{ if (NR!=1) print $3,$1,$2,$8}' AgeOfInitiation.txt | sed -e 's/ /\t/g' > temp
cat header temp > /data/neurogen/MRandPD/Results/Clumping/Data4Clump/ageofinitiation

# 3.1 COLLATE ALL SNPS FROM GWAS
cd /data/neurogen/MRandPD/Results/Clumping/Data4Clump/
for i in cigarettespday eversmoking smokingstatus ageofinitiation
do echo $i
awk '{ if (NR != 1) print $1}' $i >> allSnpRsGWASsmoking2019
done
# 3.3 SORT RS NUMBERS AND REMOVE DUPLICATED SNPS
sort allSnpRsGWASsmoking2019 | uniq > allSnpRsUniquesmoking2019
# 3.4 GET RID OF INDELS (THEY DO NOT HAVE AN RS NUMBER)
grep rs allSnpRsUniquesmoking2019 > allSnpRsUniqueNoIndelssmoking2019

# We will be ussing the 1000Genomes data in Plink format downloaded from: http://www.cog-genomics.org/plink/1.9/resources for clumping

# 4 CREATE .lsf FILES TO RUN IN CLUSTER QUEUE
#Parameters for clumping: --clump-p1 .00000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose
#this meand that intex snps will reach genomewide significance (.00000005) and either have a r2=0.001 with the next index SNP ore be 10000 kb appart
for i in cigarettespday eversmoking smokingstatus ageofinitiation
do
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
cp /data/neurogen/MRandPD/Scripts/test.lsf /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/Smoking2019ScrptClumppvalmin8Verbose 
mv /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/Smoking2019ScrptClumppvalmin8Verbose/test.lsf /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/Smoking2019ScrptClumppvalmin8Verbose/"$i".clump_chr"$j".lsf
echo "module load plink/1.90b3" >> /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/Smoking2019ScrptClumppvalmin8Verbose/"$i".clump_chr"$j".lsf
echo "plink --bfile /data/neurogen/MRandPD/1000GenomsPlinkFormat/1kg_phase1_chr"$j" --extract /data/neurogen/MRandPD/Results/Clumping/Data4Clump/allSnpRsUniqueNoIndelssmoking2019 --clump /data/neurogen/MRandPD/Results/Clumping/Data4Clump/"$i" --clump-p1 .00000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose --out /data/neurogen/MRandPD/Results/Clumping/Smoking2019Clumppvalmin8Verbose/clumped."$i".chr"$j" " >> /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/Smoking2019ScrptClumppvalmin8Verbose/"$i".clump_chr"$j".lsf
done
done

# SUBMIT THEM TO THE QUEUE
cd /data/neurogen/MRandPD/Scripts/Clmp1000GPlnkFrmt/Smoking2019ScrptClumppvalmin8Verbose/
for i in cigarettespday eversmoking smokingstatus ageofinitiation
do
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
bsub -q normal < "$i".clump_chr"$j".lsf
done 
done

#5.EXTRACT THE LIST OF INDEXSNPS FROM CLIMPING RESULTS
cd /data/neurogen/MRandPD/Results/Clumping/Smoking2019Clumppvalmin8Verbose/
cat clumped.cigarettespday.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPScigperdayverb
cat clumped.eversmoking.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPSeversmokverb
cat clumped.smokingstatus.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPSsmokingstat
cat clumped.ageofinitiation.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPSageofinit

#6.EXTRACT EXPOSURE GWAS INFO (specially BETA, MAF,SE) FOR CALCUATING PHENOTYPIC VARIANCE EXPALINED (Rsquare) AND F-STATISTIC
echo "CHROM POS RSID REF ALT AF STAT PVALUE BETA SE N EFFECTIVE_N Number_of_Studies ANNO ANNOFULL" | sed -e 's/ /\t/g'  > header
grep -w -f indexSNPSageofinit /data/neurogen/MRandPD/GWASsummaryStats/Smoking2019/AgeOfInitiation.txt | sort -n -k 2 | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticprePD/ageofinitiationforR2calc
grep -w -f indexSNPScigperdayverb /data/neurogen/MRandPD/GWASsummaryStats/Smoking2019/CigarettesPerDay.txt | sort -n -k 2 | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticprePD/cigarettespdayforR2calc
grep -w -f indexSNPSeversmokverb /data/neurogen/MRandPD/GWASsummaryStats/Smoking2019/EverSmoking.txt | sort -n -k 2 | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticprePD/eversmokingforR2calc
grep -w -f indexSNPSsmokingstat /data/neurogen/MRandPD/GWASsummaryStats/Smoking2019/SmokingStatus.txt | sort -n -k 2 | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticprePD/smokingstatusforR2calc

#7. EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET
cd /data/neurogen/MRandPD/Results/Clumping/Smoking2019Clumppvalmin8Verbose/
echo "MarkerName Allele1 Allele2 Effect StdErr P-value NStudies HetISq" | sed -e 's/ /\t/g'  > header
grep -w -f indexSNPSageofinit /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/ageofinitsmkSNPinPDdataV
grep -w -f indexSNPScigperdayverb /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/cigperdaySNPinPDdataV
grep -w -f indexSNPSeversmokverb /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/eversmkSNPinPDdataV
grep -w -f indexSNPSsmokingstat /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/smokingsatSNPinPDdataV

#8.GET MISSING VARIANTS (INDEX SNPS THAT WEREN'T FOUND ON PD GWAS DATASET)
awk '{ if (NR!=1) print $1}' ageofinitsmkSNPinPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/Clumping/Smoking2019Clumppvalmin8Verbose/indexSNPSageofinit > missingindexSNPageofinitV
awk '{ if (NR!=1) print $1}' cigperdaySNPinPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/Clumping/Smoking2019Clumppvalmin8Verbose/indexSNPScigperdayverb > missingindexSNPcigperday
awk '{ if (NR!=1) print $1}' eversmkSNPinPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/Clumping/Smoking2019Clumppvalmin8Verbose/indexSNPSeversmokverb > missingindexSNPeversmoking
awk '{ if (NR!=1) print $1}' smokingsatSNPinPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/Clumping/Smoking2019Clumppvalmin8Verbose/indexSNPSsmokingstat > missingindexSNPsmokingstat

#9.IDENTIFY PROXY SNPS FROM MISSING VARIANT (THE ONES WITH Rsquare > 0.9 WITH INDEX SNP)
cd /data/neurogen/MRandPD/Results/Clumping/Smoking2019Clumppvalmin8Verbose/
cat clumped.cigarettespday.chr*.clumped > temp
sed -e 's/     / /g' temp| sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > test
awk -F'[ ]' '$4>0.9 || $5==1' test | grep 'rs' | grep -v "^\s[0-9]" > rsqSNPscigperday
cat clumped.eversmoking.chr*.clumped > temp
sed -e 's/     / /g' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > test
awk -F'[ ]' '$4>0.9 || $5==1' test | grep 'rs' | grep -v "^\s[0-9]"> rsqSNPseversmok
cat clumped.smokingstatus.chr*.clumped > temp
sed -e 's/     / /g' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > test
awk -F'[ ]' '$4>0.9 || $5==1' test | grep 'rs' | grep -v "^\s[0-9]" > rsqSNPsmokstat
cat clumped.ageofinitiation.chr*.clumped > temp
sed -e 's/     / /g' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > test
awk -F'[ ]' '$4>0.9 || $5==1' test | grep 'rs' | grep -v "^\s[0-9]" > rsqSNPsageofinit

#10."MANUALLY" DETERMINE HOW MANY PROXY SNPS WERE FOUND DURING THE CLUMPING FOR THE MISSING SNPS
grep -w -A10 'everymissingsnp' rsqSNPs'verytrait' 

#11.EXTRACT THE PROXY SNPS THAT WERE FOUND FOR ANY OF THE MISSING SNPs
grep -A4 'rs76214862' rsqSNPseversmok | grep -v INDEX | cut -d ' ' -f2 > rs76214862RsqSNPeversmk
grep -A7 'rs12112638' rsqSNPseversmok | grep -v INDEX | cut -d ' ' -f2 > rs12112638RsqSNPeversmk

#12.EXTRAXT THIS NEW PROXY SNPs FROM THE PD GWAS DATASET, IF THERE ARE FOUND
grep -w -f rs76214862RsqSNPeversmk /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs76214862RsqSNPeversmkinpPD
grep -w -f rs12112638RsqSNPeversmk /data/neurogen/MRandPD/GWASsummaryStats/pd/META_ANALYSIS_10K23_beta_se_wo23AndMe1_pdgene_sharing_tab.tbl | sort | uniq > /data/neurogen/MRandPD/Results/HarmonizationVerbose/rs12112638RsqSNPeversmkinpPD

#13.JOIN ALL INDEX SNPS FOUND ON THE PD GWAS DATASET
cat eversmkSNPinPDdataV rs76214862RsqSNPeversmkinpPD rs12112638RsqSNPeversmkinpPD > fulleversmkSNPinPDdataV

#EXTRACT THE INFO FOR CALCULATING Rsquare AND F-STATISTIC FROM INDEX SNPS ALSO PRESENT IN THE PD GWAS DATASET 
cd /data/neurogen/MRandPD/Results/HarmonizationVerbose/
echo "CHROM POS RSID REF ALT AF STAT PVALUE BETA SE N EFFECTIVE_N Number_of_Studies ANNO ANNOFULL" | sed -e 's/ /\t/g'  > header
awk '{ if (NR!=1) print $1}' ageofinitsmkSNPinPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/Smoking2019/AgeOfInitiation.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/ageofinitiationforR2calcV
awk '{ if (NR!=1) print $1}' cigperdaySNPinPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/Smoking2019/CigarettesPerDay.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/cigarettespdayforR2calcV
awk '{ if (NR!=1) print $1}' fulleversmkSNPinPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/Smoking2019/EverSmoking.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/eversmokingforR2calcV
awk '{ if (NR!=1) print $1}' smokingsatSNPinPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/Smoking2019/SmokingStatus.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/smokingstatusforR2calcV
