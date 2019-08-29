##This code is for obtainig the index SNPs (IV) by clumping from smoking GWAS data
##We need to clump SNPs to ensure that all the index SNPs are independent frome one another
##After that, we will extract those index SNPs data from the PD GWAS, an calculate the Rsquare an F-statistic from the instruments
##We are using the summary statistics from Liu,2019 downloaded from https://conservancy.umn.edu/handle/11299/201564

#2. PREPARE GWAS FILES FOR SNP CLUMPING (in PLINK)
# EXTRACTION OF COLUMNS (SNP CHR BP P) from the GWAS results
mkdir /data/neurogen/MRandPD/Results/PLINKClumping
mkdir /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping
cd /data/neurogen/MRandPD/GWASsummaryStats/Smoking2019/
echo "SNP CHR BP P" | sed -e 's/ /\t/g'  > header
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation
do 
awk '{ if (NR!=1) print $3,$1,$2,$8}' $i.txt | sed -e 's/ /\t/g' > temp
cat header temp > /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/$i
done

# 3.1 COLLATE ALL SNPS FROM GWAS
cd /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation
do echo $i
awk '{ if (NR != 1) print $1}' $i >> allSnpRsGWASSmoking2019
done
# 3.3 SORT RS NUMBERS AND REMOVE DUPLICATED SNPS
sort allSnpRsGWASSmoking2019 | uniq > allSnpRsUniqueSmoking2019
# 3.4 GET RID OF INDELS (THEY DO NOT HAVE AN RS NUMBER)
grep rs allSnpRsUniqueSmoking2019 > allSnpRsUniqueNoIndelsSmoking2019

# We will be ussing the 1000Genomes data in Plink format downloaded from: http://www.cog-genomics.org/plink/1.9/resources for clumping

# 4 CREATE .lsf FILES TO RUN IN CLUSTER QUEUE
#Parameters for clumping: --clump-p1 .00000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose
#this mean that index snps will reach genomewide significance (.00000005) and either have a r2=0.001 with the next index SNP or be 10000 kb appart
#--clump-verbose to ger the Rsquare values of the clumped SNPs with the index SNP
mkdir /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/
mkdir /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Smoking2019ScrptClumpPmin8V/
mkdir /data/neurogen/MRandPD/Results/PLINKClumping/Smoking2019ClumpPmin8V/
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation
do
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
cp /data/neurogen/MRandPD/Scripts/test.lsf /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Smoking2019ScrptClumpPmin8V 
mv /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Smoking2019ScrptClumpPmin8V/test.lsf /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Smoking2019ScrptClumpPmin8V/"$i".clump_chr"$j".lsf
echo "module load plink/1.90b3" >> /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Smoking2019ScrptClumpPmin8V/"$i".clump_chr"$j".lsf
echo "plink --bfile /data/neurogen/MRandPD/1000GenomsPlinkFormat/1kg_phase1_chr"$j" --extract /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/allSnpRsUniqueNoIndelsSmoking2019 --clump /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/"$i" --clump-p1 .00000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose --out /data/neurogen/MRandPD/Results/PLINKClumping/Smoking2019ClumpPmin8V/clumped."$i".chr"$j" " >> /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Smoking2019ScrptClumpPmin8V/"$i".clump_chr"$j".lsf
done
done

# SUBMIT THEM TO THE QUEUE
cd /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/Smoking2019ScrptClumpPmin8V/
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation
do
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
bsub -q normal < "$i".clump_chr"$j".lsf
done 
done

#Verify that all clumping was compleated. to do so, run the following command and it whould return 88. if not, look for the ones missing and rerun them 
cd /data/neurogen/MRandPD/Results/PLINKClumping/Smoking2019ClumpPmin8V/
ls | grep 'log' | wc -l

#5.EXTRACT THE LIST OF INDEX SNPS FROM CLUMPING RESULTS
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation
do
cat clumped.$i.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPs"$i"V
done

#6.EXTRACT EXPOSURE GWAS INFO FOR CALCUATING PHENOTYPIC VARIANCE EXPALINED (Rsquare) AND F-STATISTIC
mkdir /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/
echo "CHROM POS RSID REF ALT AF STAT PVALUE BETA SE N EFFECTIVE_N Number_of_Studies ANNO ANNOFULL" | sed -e 's/ /\t/g'  > header
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation
do
grep -w -f indexSNPs"$i"V /data/neurogen/MRandPD/GWASsummaryStats/Smoking2019/$i.txt | sort -n -k 2 | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/"$i"forR2calc
done

#7. AS PD DATASET DO NOT HAVE rsID BUT CHR:BP VALUES, EXTRACT THE CHR:POS OF OUR INDEX SNPS
cd /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/
echo "SNP rsID" > header2
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation
do
awk '{ if (NR!=1) print $3}' "$i"forR2calc > "$i"rsID
awk '{ if (NR!=1) print $1,$2}' "$i"forR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > "$i"chrpos
paste "$i"chrpos "$i"rsID | sort > temp2
cat header2 temp2 | sed -e 's/ /\t/g' > "$i"chrposrsID
done

#8. EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET
mkdir /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/
echo "SNP A1 A2 freq b se p N_cases N_controls" > header
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation
do
grep -w -f "$i"chrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab| sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/"$i"SNPsinNallsPDdataV
done 

#9. ADD THE rsID NUMBER TO THE PD DATA
cd /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation
do
join "$i"SNPsinNallsPDdataV /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/"$i"chrposrsID | sed -e 's/ /\t/g' > temp
mv temp "$i"SNPsinNallsPDdataV
done 
