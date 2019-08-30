##EDUCATIONAL ATTAINMENT

#2. PREPARE GWAS FILES FOR SNP CLUMPING (in PLINK)
# EXTRACTION OF COLUMNS (SNP EA NEA P BETA) from the GWAS results
# eduyears
cd /data/neurogen/MRandPD/GWASsummaryStats/education/
echo "SNP CHR BP P" | sed -e 's/ /\t/g'  > header
awk '{ if (NR!=1) print $1,$2,$3,$9}' EduYears_Main.txt | sed -e 's/ /\t/g' > temp
cat header temp > /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/eduyears

# 3.1 COLLATE ALL SNPS FROM GWAS
cd /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/
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
mkdir /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/EduyearsScrptClumpPmin8V/
mkdir /data/neurogen/MRandPD/Results/PLINKClumping/EduyearsClumpPmin8V/
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
cp /data/neurogen/MRandPD/Scripts/test.lsf /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/EduyearsScrptClumpPmin8V 
mv /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/EduyearsScrptClumpPmin8V/test.lsf /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/EduyearsScrptClumpPmin8V/eduyears.clump_chr"$j".lsf
echo "module load plink/1.90b3" >> /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/EduyearsScrptClumpPmin8V/eduyears.clump_chr"$j".lsf
echo "plink --bfile /data/neurogen/MRandPD/1000GenomsPlinkFormat/1kg_phase1_chr"$j" --extract /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/allSnpRsUniqueNoIndelseduyears --clump /data/neurogen/MRandPD/Results/PLINKClumping/GWASData4Clumping/eduyears --clump-p1 .00000005 --clump-p2 .01 --clump-r2 0.001 --clump-kb 10000 --clump-verbose --out /data/neurogen/MRandPD/Results/PLINKClumping/EduyearsClumpPmin8V/clumped.eduyears.chr"$j" " >> /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/EduyearsScrptClumpPmin8V/eduyears.clump_chr"$j".lsf
done

# SUBMIT THEM TO THE QUEUE
cd /data/neurogen/MRandPD/Scripts/PLINKClumpScripts/EduyearsScrptClumpPmin8V/
for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
bsub -q normal < eduyears.clump_chr"$j".lsf
done 

#Verify that all clumping was compleated. to do so, run the following command and it whould return 22. if not, look for the ones missing and rerun them 
cd /data/neurogen/MRandPD/Results/PLINKClumping/EduyearsClumpPmin8V/
ls | grep 'log' | wc -l

#5.EXTRACT THE LIST OF INDEXSNPS FROM CLUMPING RESULTS
cat clumped.eduyears.chr*.clumped > temp
grep 'INDEX' temp | sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | cut -d ' ' -f3 > indexSNPSeduyearsV

#6.EXTRACT EXPOSURE GWAS INFO (specially BETA, MAF,SE) FOR CALCUATING PHENOTYPIC VARIANCE EXPALINED (Rsquare) AND F-STATISTIC
echo "MarkerName CHR POS A1 A2 EAF Beta SE Pval" | sed -e 's/ /\t/g'  > header
grep -w -f indexSNPSeduyearsV /data/neurogen/MRandPD/GWASsummaryStats/education/EduYears_Main.txt | sort -n -k 2 | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/eduyearsforR2calc

#7. AS PD NALLS DATASET DOESN'T HAVE AN rsID, EXTRACT THE CHR:BP POSITION
cd /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/
awk '{ if (NR!=1) print $1}' eduyearsforR2calc > eduyearsrsID
awk '{ if (NR!=1) print $2,$3}' eduyearsforR2calc > temp
sed -e 's/^/chr/g' temp | sed -e 's/ /:/g' > eduyearschrpos
paste eduyearschrpos eduyearsrsID | sort > temp2
cat header2 temp2 | sed -e 's/ /\t/g' > eduyearschrposrsID

#8. EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET
echo "SNP A1 A2 freq b se p N_cases N_controls" > header
grep -w -f eduyearschrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort | uniq > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/eduyearsSNPinNallsPDdataV

#9. GET MISSING VARIANTS (INDEX SNPS THAT WEREN'T FOUND ON PD GWAS DATASET)
cd /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/
awk '{ if (NR!=1) print $1}' eduyearsSNPinNallsPDdataV > temp
grep -v -f temp /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/eduyearschrpos | sort > missingeduyearsSNPinNallsPDchrpos
grep -f missingeduyearsSNPinNallsPDchrpos /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/eduyearschrposrsID> missingeduyearsSNPinNallsPDchrposrsID
cut -f2 missingeduyearsSNPinNallsPDchrposrsID > missingeduyearsSNPinNallsPDrsID

#10.IDENTIFY PROXY SNPS FROM MISSING VARIANT (THE ONES WITH Rsquare > 0.9 WITH INDEX SNP)
cd /data/neurogen/MRandPD/Results/PLINKClumping/EduyearsClumpPmin8V/
cat clumped.eduyears.chr*.clumped > temp
sed -e 's/     / /g' temp| sed -e 's/    / /g' | sed -e 's/   / /g' | sed -e 's/  / /g' | grep -v "dataset" > prueba
awk -F'[ ]' '$4>0.6 || $5==1' prueba | grep 'rs' | grep -v "^\s[0-9]" > rsq6SNPseduyears

#extracto all the proxys with Rsquare > 0.6
missingeduyears=/data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/missingeduyearsSNPinNallsPDrsID
rsIDmissedu=(`awk '{print}' $missingeduyears`)
for i in "${rsIDmissedu[@]}"; 
do 
sed -e "s/(INDEX)\s$i/$i/" rsq6SNPseduyears > temp 
sed -n "/$i/,/INDEX/p" temp | grep -v "INDEX" | grep -v "$i" | sort -k1 | sed -e 's/^ //g' | sed -e 's/ /\t/g' > proxysof"$i"_eduyears
cut -f1 proxysof"$i"_eduyears > rsID"$i"_eduyears
done

wc -l rsID*

cut -f1,2,3 /data/neurogen/MRandPD/GWASsummaryStats/education/EduYears_Main.txt > eduyearschrposID

#once we have extracted them, we have to convert theri rsID to ch:pos so we can look for them on the PD dataset
for i in "${rsIDmissedu[@]}"
do
grep -w -f rsID"$i"_eduyears eduyearschrposID > temp 
cut -f1 temp > rsid
cut -f2,3 temp | sed -e 's/^/chr/g' | sed -e 's/\t/:/g' > chrpos
paste rsid chrpos | sort -k2 > "$i"chrposrsID
cut -f2 "$i"chrposrsID > "$i"chrpos
done

#11. EXTRACT THE INFORMATION OF THE INDEX SNPS BUT FROM THE PD GWAS DATASET
for i in "${rsIDmissedu[@]}"
do
grep -w -f "$i"chrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort -k1 > temp
cut -f1 temp > temp2 
grep -w -f temp2 "$i"chrposrsID | sort -k2 > temp3
join -1 1 -2 2 temp temp3 | sed -e 's/ /\t/g' | sort -k10 > temp4
cut -f10 temp4 > temp5 
grep -w -f temp5 proxysof"$i"_eduyears | sort -k1 > temp6 
join -1 10 -2 1 temp4 temp6 | sort -k15,15g -k12,12nr | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$1}' | head -n1 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/"$i"eduyears_SNPsinNallsPDdataV
done

grep -w -f "$i"chrpos /data/neurogen/MRandPD/GWASsummaryStats/PDnalls2019/nallsEtAl2019_excluding23andMe_allVariants.tab | sort -k1 > temp
cut -f1 temp > temp2 
grep -w -f temp2 "$i"chrposrsID | sort -k2 > temp3
join -1 1 -2 2 temp temp3 | sed -e 's/ /\t/g' | sort -k10 > temp4
cut -f10 temp4 > temp5 
grep -w -f temp5 proxysof"$i"_eduyears | sort -k1 > temp6 
join -1 10 -2 1 temp4 temp6 | sort -k15,15g -k12,12nr | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$1}' | head -n1 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/"$i"eduyears_SNPsinNallsPDdataV


cd /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/
wc -l *eduyears_SNPsinNallsPDdataV

#12. ADD THE rsID NUMBER TO THE PD DATA
cd /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/
join eduyearsSNPinNallsPDdataV /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/eduyearschrposrsID | sed -e 's/ /\t/g' > temp
cat *eduyears_SNPsinNallsPDdataV  > temp2
cat temp temp2 | sed -e 's/ /\t/g' > eduyearsSNPinNallsPDdataV

#13. THEN WE WILL EXTRACT THE INFORMATION FOR THIS NEW SET TO SNPs TO CALCULATE THEIR Rsq AND F-STATISTIC
echo "SNP CHR BP A1 A2 MAF BETA SE P N" | sed -e 's/ /\t/g'  > header
awk '{ if (NR!=1) print $10}' eduyearsSNPinNallsPDdataV > temp
grep -w -f temp /data/neurogen/MRandPD/GWASsummaryStats/education/EduYears_Main.txt | sort -n -k 2 | uniq > temp2
cat header temp2 | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/eduyearsforR2calcVNalls
