##This code is for obtainig the index SNPs (IV) by clumping from SLEEPING GWAS data
##We need to clump SNPs to ensure that all the index SNPs are independent frome one another
##After that, we will extract those index SNPs data from the PD GWAS, an calculate the Rsquare an F-statistic from the instruments
##We are using the summary statistics from Lane,2017 https://doi.org/10.1038/ng.3749

# 2. PREPARE GWAS FILES FOR SNP CLUMPING (in PLINK)
# EXTRACTION OF COLUMNS (SNP EA NEA P BETA) from the GWAS results
#Issue with chromosome 23? in sleeping disturbances data
#get rid of the extra info following the snp rs number, select just snp bp pval columns 
cd /data/neurogen/MRandPD/GWASsummaryStats/sleep/
mv Application6818_SleepTrait_data_Lane_RESULTS_DATA_EXCESSIVE_DAYTIME_SLEEPINESS_all.txt Application6818_SleepTrait_data_Lane_RESULTS_DATA_EXCESSDAYSLEEP_all.txt
mv Application6818_SleepTrait_data_Lane_RESULTS_DATA_INSOMNIA_SYMPTOMS_all.txt Application6818_SleepTrait_data_Lane_RESULTS_DATA_INSOMNIA_all.txt
mv Application6818_SleepTrait_data_Lane_RESULTS_DATA_SLEEP_DURATION_all.txt Application6818_SleepTrait_data_Lane_RESULTS_DATA_SLEEPDUR_all.txt
echo "SNP BP P" | sed -e 's/ /\t/g'  > header
for i in EXCESSDAYSLEEP INSOMNIA SLEEPDUR 
do 
awk '{ if (NR!=1) print $1,$3,$9}' Application6818_SleepTrait_data_Lane_RESULTS_DATA_"$i"_all.txt | sed -e 's/ /\t/g' > temp
grep ':' temp  > issue
grep -F -x -v -f issue temp > data
sed -e 's/:/\t/g' issue | cut -f1,5,6 > temp
cat header data temp > /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/"$i"wochr23
done

# 3.1 COLLATE ALL SNPS FROM GWAS
# For sleeping GWAS
cd /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/
for i in EXCESSDAYSLEEP INSOMNIA SLEEPDUR 
do echo $i
awk '{ if (NR != 1) print $1}' "$i"wochr23 >> allSnpGWASSleepwochr23
done
# 3.3 SORT RS NUMBERS AND REMOVE DUPLICATED SNPS
sort allSnpGWASSleepwochr23 | uniq > allSnpGWASUniqueSleepwochr23
# 3.4 GET RID OF INDELS (THEY DO NOT HAVE AN RS NUMBER)
grep 'rs' allSnpGWASUniqueSleepwochr23 > allSnpGWASUniqueNoIndelsSleepwochr23

# We will be ussing the 1000Genomes data in Plink format downloaded from: http://www.cog-genomics.org/plink/1.9/resources for clumping

# 4 CREATE .lsf FILES TO RUN IN CLUSTER QUEUE
#Parameters for clumping: --clump-p1 .0000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose
#for identifying both genomewide significant and suggestive index snps which are independent form one another
#this mean that intex snps will reach significance (.0000005) and either have a r2=0.001 with the next index SNP or be 10000 kb appart
#--clump-verbose to ger the Rsquare values of the clumped SNPs with the index SNP
mkdir /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/SleepScrptClumpPmin7V/
mkdir /data/neurogen/MRandPD/Results/PLINKClumping/SleepClumpPmin7V/
for i in EXCESSDAYSLEEP INSOMNIA SLEEPDUR
do
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
cp /data/neurogen/MRandPD/Scripts/test.lsf /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/SleepScrptClumpPmin7V/ 
mv /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/SleepScrptClumpPmin7V/test.lsf /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/SleepScrptClumpPmin7V/"$i".clump_chr"$j".lsf
echo "module load plink/1.90b3" >> /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/SleepScrptClumpPmin7V/"$i".clump_chr"$j".lsf
echo "plink --bfile /data/neurogen/MRandPD/1000GenomsPlinkFormat/1kg_phase1_chr"$j" --extract /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/allSnpGWASUniqueNoIndelsSleepwochr23 --clump /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/"$i"wochr23 --clump-p1 .0000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose --out /data/neurogen/MRandPD/Results/PLINKClumping/SleepClumpPmin7V/clumped."$i".chr"$j" " >> /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/SleepScrptClumpPmin7V/"$i".clump_chr"$j".lsf
done
done

# SUBMIT THEM TO THE QUEUE
cd /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/SleepScrptClumpPmin7V/
for i in  EXCESSDAYSLEEP INSOMNIA SLEEPDUR
do
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 
do
bsub -q normal < "$i".clump_chr"$j".lsf
done
done

#Verify that all clumping was compleated. to do so, run the following command and it whould return 66. if not, look for the ones missing and rerun them 
cd /data/neurogen/MRandPD/Results/PLINKClumping/SleepClumpPmin7V/
ls | grep 'log' | wc -l

#5.EXTRACT THE LIST OF INDEXSNPS FROM CLUMPING RESULTS
for i in  EXCESSDAYSLEEP INSOMNIA SLEEPDUR
do
cat clumped.$i.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPs"$i"V
done

#6.EXTRACT EXPOSURE GWAS INFO (specially BETA, MAF,SE) FOR CALCUATING PHENOTYPIC VARIANCE EXPALINED (Rsquare) AND F-STATISTIC
echo "SNP CHR BP A1 A2 MAF BETA SE P N" | sed -e 's/ /\t/g'  > header
for i in  EXCESSDAYSLEEP INSOMNIA SLEEPDUR
do
grep -w -f indexSNPs"$i"V /data/neurogen/MRandPD/GWASsummaryStats/sleep/Application6818_SleepTrait_data_Lane_RESULTS_DATA_"$i"_all.txt | sort -n -k 2 | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/"$i"forR2calc
done

#7. AS PD DATASET DO NOT HAVE rsID BUT CHR:BP VALUES, EXTRACT THE CHR:POS OF OUR INDEX SNPS
cd /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/
for i in EXCESSDAYSLEEP INSOMNIA SLEEPDUR
do
awk '{ if (NR!=1) print $1}' "$i"forR2calc > "$i"rsID
awk '{ if (NR!=1) print $2,$3}' "$i"forR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > "$i"chrpos
paste "$i"chrpos "$i"rsID | sort > temp2
cat header2 temp2 | sed -e 's/ /\t/g' > "$i"chrposrsID
done

#8. EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET
for i in EXCESSDAYSLEEP INSOMNIA SLEEPDUR
do
grep -w -f "$i"chrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/"$i"SNPsinNallsPDdataV
done 

#10.GET MISSING VARIANTS (INDEX SNPS THAT WEREN'T FOUND ON PD GWAS DATASET)
cd /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/
for i in EXCESSDAYSLEEP INSOMNIA SLEEPDUR
do
awk '{ if (NR!=1) print $1}' "$i"SNPsinNallsPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/"$i"chrpos | sort > missing"$i"SNPinNallsPDchrpos
grep -f missing"$i"SNPinNallsPDchrpos /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/"$i"chrposrsID > missing"$i"SNPinNallsPDchrposrsID
cut -f2 missing"$i"SNPinNallsPDchrposrsID > missing"$i"SNPinNallsPDrsID
done 

#11.IDENTIFY PROXY SNPS FROM MISSING VARIANT (THE ONES WITH Rsquare > 0.6 WITH INDEX SNP)
cd /data/neurogen/MRandPD/Results/PLINKClumping/SleepClumpPmin7V/
for i in EXCESSDAYSLEEP INSOMNIA SLEEPDUR
do
cat clumped.$i.chr*.clumped > temp
sed -e 's/     / /g' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > test
awk -F'[ ]' '$4>0.6 || $5==1' test | grep 'rs' | grep -v "^\s[0-9]" | sed -e 's/^ //g' |  sed -e 's/ /\t/g' > rsq6SNPs"$i"
done 

sed -e "s/(INDEX)\s$i/$i/" rsq6SNPsINSOMNIA > temp 
sed -n "/$i/,/INDEX/p" temp | grep -v "INDEX" | grep -v "$i" | sort -k1 > proxysof"$i"_INSOMNIA
#no proxy was found for the missing insomnia SNP

#extracto all the proxys with Rsquare > 0.6
missingEXCESSDAYSLEEP=/data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/missingEXCESSDAYSLEEPSNPinNallsPDrsID
rsIDmissEDS=(`awk '{print}' $missingEXCESSDAYSLEEP`)
for i in "${rsIDmissEDS[@]}"; 
do 
sed -e "s/(INDEX)\s$i/$i/" rsq6SNPsEXCESSDAYSLEEP > temp 
sed -n "/$i/,/INDEX/p" temp | grep -v "INDEX" | grep -v "$i" | sort -k1 > proxysof"$i"_EXCESSDAYSLEEP
cut -f1 proxysof"$i"_EXCESSDAYSLEEP > rsID"$i"_EXCESSDAYSLEEP
done

wc -l rsID*_EXCESSDAYSLEEP
#proxy snps were found for 5 snps

missingSLEEPDUR=/data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/missingSLEEPDURSNPinNallsPDrsID
rsIDmissSDur=(`awk '{print}' $missingSLEEPDUR`)
for i in "${rsIDmissSDur[@]}"; 
do 
sed -e "s/(INDEX)\s$i/$i/" rsq6SNPsSLEEPDUR > temp 
sed -n "/$i/,/INDEX/p" temp | grep -v "INDEX" | grep -v "$i" > proxysof"$i"_SLEEPDUR
cut -f1 proxysof"$i"_SLEEPDUR > rsID"$i"_SLEEPDUR
done

wc -l proxy*_SLEEPDUR
#no proxy snps were found for any of the missing snps

find . -size 0 -delete

cut -f1,2,3  -d ' ' /data/neurogen/MRandPD/GWASsummaryStats/sleep/Application6818_SleepTrait_data_Lane_RESULTS_DATA_EXCESSDAYSLEEP_all.txt > ESCESSDAYSLEEPchrposID

ls rsID*_EXCESSDAYSLEEP | sed -e 's/rsID//g' | sed -e 's/_EXCESSDAYSLEEP//g' > allproxysEDS
allEDSproxys=allproxysEDS
EDSallproxy=(`awk '{print}' $allEDSproxys`)

#once we have extracted them, we have to convert theri rsID to ch:pos so we can look for them on the PD dataset
for i in "${EDSallproxy[@]}"
do
grep -w -f rsID"$i"_EXCESSDAYSLEEP ESCESSDAYSLEEPchrposID > temp 
cut -f1 -d ' ' temp > rsid
cut -f2,3 -d ' ' temp | sed -e 's/^/chr/g' | sed -e 's/ /:/g' > chrpos
paste rsid chrpos | sort -k2 > "$i"chrposrsID
cut -f2 "$i"chrposrsID > "$i"chrpos
done

#8. EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET

#now we are extracting the proxy that was found on the PD dataset and has the best scores for both the pvalue and the rsquare
for i in "${EDSallproxy[@]}"
do
grep -w -f "$i"chrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort -k1 > temp
cut -f1 temp > temp2 
grep -w -f temp2 "$i"chrposrsID | sort -k2 > temp3
join -1 1 -2 2 temp temp3 | sed -e 's/ /\t/g' | sort -k10 > temp4
cut -f10 temp4 > temp5 
grep -w -f temp5 proxysof"$i"_EXCESSDAYSLEEP | sort -k1 > temp6 
join -1 10 -2 1 temp4 temp6 | sort -k15,15n -k13,13nr | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$1}' | head -n1 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/"$i"_SNPsinNallsPDdataV
done

cd /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/
find . -size 0 -delete

#9. ADD THE rsID NUMBER TO THE PD DATA
cd /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/
for i in EXCESSDAYSLEEP INSOMNIA SLEEPDUR
do
join "$i"SNPsinNallsPDdataV /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/"$i"chrposrsID | sed -e 's/ /\t/g' > temp
mv temp "$i"SNPsinNallsPDdataV
done

cat EXCESSDAYSLEEPSNPsinNallsPDdataV rs115320831_SNPsinNallsPDdataV rs35309287_SNPsinNallsPDdataV > temp
mv temp EXCESSDAYSLEEPSNPsinNallsPDdataV

#16. THEN WE WILL EXTRACT THE INFORMATION FOR THIS NEW SET TO SNPs TO CALCULATE THEIR Rsq AND F-STATISTIC
echo "SNP CHR BP A1 A2 MAF BETA SE P N" | sed -e 's/ /\t/g'  > header
for i in EXCESSDAYSLEEP INSOMNIA SLEEPDUR
do
awk '{ if (NR!=1) print $10}' "$i"SNPsinNallsPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/sleep/Application6818_SleepTrait_data_Lane_RESULTS_DATA_EXCESSDAYSLEEP_all.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/"$i"forR2calcVNalls
done
