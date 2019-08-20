##This code is for obtainig the index SNPs (IV) by clumping from smoking GWAS data
##We need to clump SNPs to ensure that all the index SNPs are independent frome one another
##After that, we will extract those index SNPs data from the PD GWAS, an calculate the Rsquare an F-statistic from the instruments
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
#this mean that index snps will reach genomewide significance (.00000005) and either have a r2=0.001 with the next index SNP or be 10000 kb appart
#--clump-verbose to ger the Rsquare values of the clumped SNPs with the index SNP
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

#5.EXTRACT THE LIST OF INDEXSNPS FROM CLUMPING RESULTS
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

#7. AS PD DATASET DO NOT HAVE rsID BUT CHR:BP VALUES, EXTRACT THE CHR:POS OF OUR INDEX SNPS
awk '{ if (NR!=1) print $1,$2}' ageofinitiationforR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > ageofinitchrpos
awk '{ if (NR!=1) print $1,$2}' cigarettespdayforR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > cigarettespdaychrpos
awk '{ if (NR!=1) print $1,$2}' eversmokingforR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > eversmokingchrpos
awk '{ if (NR!=1) print $1,$2}' smokingstatusforR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > smokingstatchrpos

awk '{ if (NR!=1) print $3}' cigarettespdayforR2calc > cigarettespdayrsID
paste cigarettespdaychrpos cigarettespdayrsID | sort > temp
cat header2 temp | sed -e 's/ /\t/g' > cigarettespdaychrposrsID
awk '{ if (NR!=1) print $3}' drinksperweekforR2calc > drinkspweekrsID
paste drinkspweekchrpos drinkspweekrsID | sort > temp
cat header2 temp | sed -e 's/ /\t/g' > drinkspweekchrposrsID
awk '{ if (NR!=1) print $3}' eversmokingforR2calc > eversmokingrsID
paste eversmokingchrpos eversmokingrsID | sort > temp
cat header2 temp | sed -e 's/ /\t/g' > eversmokingposrsID
awk '{ if (NR!=1) print $3}' smokingstatusforR2calc > smokingstatrsID
paste smokingstatchrpos smokingstatrsID | sort > temp
cat header2 temp | sed -e 's/ /\t/g' > smokingstatposrsID

#7. EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET
echo "SNP A1 A2 freq b se p N_cases N_controls" > header
grep -w -f ageofinitchrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab| sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/ageofinitSNPinNallsPDdataV
grep -w -f cigarettespdaychrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/cigarettespdaySNPinNallsPDdataV
grep -w -f eversmokingchrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab| sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/eversmokingSNPinNallsPDdataV
grep -w -f smokingstatchrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/smokingstatSNPinNallsPDdataV

#AD THE rsID NUMBER TO THE DATA
join ageofinitSNPinNallsPDdataV /data/neurogen/MRandPD/Results/R2andFstatisticprePD/ageofinitchrposrsID | sed -e 's/ /\t/g' > temp
mv temp ageofinitSNPinNallsPDdataV
join cigarettespdaySNPinNallsPDdataV /data/neurogen/MRandPD/Results/R2andFstatisticprePD/cigarettespdaychrposrsID | sed -e 's/ /\t/g' > temp
mv temp cigarettespdaySNPinNallsPDdataV
join drinkspweekSNPinNallsPDdataV /data/neurogen/MRandPD/Results/R2andFstatisticprePD/drinkspweekchrposrsID | sed -e 's/ /\t/g' > temp
mv temp drinkspweekSNPinNallsPDdataV
join eversmokingSNPinNallsPDdataV /data/neurogen/MRandPD/Results/R2andFstatisticprePD/eversmokingposrsID | sed -e 's/ /\t/g' > temp
mv temp eversmokingSNPinNallsPDdataV
join smokingstatSNPinNallsPDdataV /data/neurogen/MRandPD/Results/R2andFstatisticprePD/smokingstatposrsID | sed -e 's/ /\t/g' > temp
mv temp smokingstatSNPinNallsPDdataV
