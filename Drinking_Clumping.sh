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
cat header temp > /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/drinksperweek

# 3.1 COLLATE ALL SNPS FROM GWAS
cd /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/
awk '{ if (NR != 1) print $1}' drinksperweek > allSnpRsGWASdrinking
# 3.3 SORT RS NUMBERS AND REMOVE DUPLICATED SNPS
sort allSnpRsGWASdrinking | uniq > allSnpRsUniquedrinking
# 3.4 GET RID OF INDELS (THEY DO NOT HAVE AN RS NUMBER)
grep rs allSnpRsUniquedrinking > allSnpRsUniqueNoIndelsdrinking

# Using the 1000Genomes data in Plink format downloaded from: http://www.cog-genomics.org/plink/1.9/resources

# 4. CREATE .lsf FILES TO RUN IN CLUSTER QUEUE
#Parameters for clumping: --clump-p1 .00000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose
#this mean that index snps will reach genomewide significance (.00000005) and either have a r2=0.001 with the next index SNP or be 10000 kb appart
#--clump-verbose to ger the Rsquare values of the clumped SNPs with the index SNP
mkdir /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/DrinkingScrptClumpPmin8V/
mkdir /data/neurogen/MRandPD/Results/PLINKClumping/DrinkingClumpPmin8V/
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
cp /data/neurogen/MRandPD/Scripts/test.lsf /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/DrinkingScrptClumpPmin8V/ 
mv /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/DrinkingScrptClumpPmin8V/test.lsf /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/DrinkingScrptClumpPmin8V/drinksperweek.clump_chr"$j".lsf
echo "module load plink/1.90b3" >> /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/DrinkingScrptClumpPmin8V/drinksperweek.clump_chr"$j".lsf
echo "plink --bfile /data/neurogen/MRandPD/1000GenomsPlinkFormat/1kg_phase1_chr"$j" --extract /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/allSnpRsUniqueNoIndelsdrinking --clump /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/drinksperweek --clump-p1 .00000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose --out /data/neurogen/MRandPD/Results/PLINKClumping/DrinkingClumpPmin8V/clumped.drinksperweek.chr"$j" " >> /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/DrinkingScrptClumpPmin8V/drinksperweek.clump_chr"$j".lsf
done

# SUBMIT THEM TO THE QUEUE
cd /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/DrinkingScrptClumpPmin8V/
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
bsub -q normal < drinksperweek.clump_chr"$j".lsf
done 

#Verify that all clumping was compleated. to do so, run the following command and it whould return 22. if not, look for the ones missing and rerun them 
cd /data/neurogen/MRandPD/Results/PLINKClumping/DrinkingClumpPmin8V/
ls | grep 'log' | wc -l

#5.EXTRACT THE LIST OF INDEXSNPS FROM CLUMPING RESULTS
cat clumped.drinksperweek.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPSdrinksperweekV

#6.EXTRACT EXPOSURE GWAS INFO (specially BETA, MAF,SE) FOR CALCUATING PHENOTYPIC VARIANCE EXPALINED (Rsquare) AND F-STATISTIC
echo "CHROM POS RSID REF ALT AF STAT PVALUE BETA SE N EFFECTIVE_N Number_of_Studies ANNO ANNOFULL" | sed -e 's/ /\t/g'  > header
grep -w -f indexSNPSdrinksperweekV /data/neurogen/MRandPD/GWASsummaryStats/drinking/DrinksPerWeek.txt | sort -n -k 2 | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/drinksperweekforR2calc

#7. AS PD NALLS DATASET DOESN'T HAVE AN rsID, EXTRACT THE CHR:BP POSITION
cd /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/
awk '{ if (NR!=1) print $1,$2}' drinksperweekforR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > drinkspweekchrpos
awk '{ if (NR!=1) print $3}' drinksperweekforR2calc > drinkspweekrsID
paste drinkspweekchrpos drinkspweekrsID | sort > temp
cat header2 temp | sed -e 's/ /\t/g' > drinkspweekchrposrsID

#8. EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET
echo "SNP A1 A2 freq b se p N_cases N_controls" > header
grep -w -f drinkspweekchrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/drinkspweekSNPinNallsPDdataV

#9. AD THE rsID
cd /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/
join drinkspweekSNPinNallsPDdataV /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/drinkspweekchrposrsID | sed -e 's/ /\t/g' > temp
mv temp drinkspweekSNPinNallsPDdataV
