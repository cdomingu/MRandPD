#This Script is for conducting Harmonization to ensure that both datasets that contain the indexSNPs for each exposure (IV-exposure, IV-outcome),
#are identically coded regarding the effect allele. If not, we will flip the one from the PD GWAS Dataset so that it matches its exposure counterpart


#2. Extract the information from the SNPs present in the PD GWAS data for each exposure
cd /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/
mkdir /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/
mkdir /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/InstruVar/
mkdir /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/IVinHarmonyPositive/
echo "RsID REF ALT AF BETA SE PVAL N" > header
for i in EXCESSDAYSLEEP INSOMNIA SLEEPDUR eduyears
do
awk '{ if (NR!=1) print $1,$4,$5,$6,$7,$8,$9,$10}' /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/AfterNallsR2calc/"$i"forR2calcVNalls > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/InstruVar/"$i"IV
done

for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation drinksperweek
do
awk '{ if (NR!=1) print $3,$4,$5,$6,$9,$10,$8,$11}' /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/smoking4R2calc/"$i"forR2calc > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/InstruVar/"$i"IV
done
#all the files in the directory /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/InstruVar/
# will be the input for the first step of the Rscript titled Harmonization.Rmd, which will returned the *fliped files

echo "RsID REF OR PVAL N" > header2
for i in blonde brown red
do
awk '{ if (NR!=1) print $2,$4,$7,$9,$6}' /data/neurogen/MRandPD/Results/R2andFstat4IndexSNPs/hair4R2calc/"$i"forR2calc > temp
cat header2 temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/IVinHarmonyPositive/"$i"
done

#3. AFTER HARMONIZING THE IV TO ALL POSITIVE EFFECT SIZE (Harmonizarion.rmd) AND EXTRACT THE EFFECT ALLLELE CODING
cd /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/IVinHarmonyPositive/
mkdir /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/IVEA/
echo "RsID REF" > header
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation drinksperweek EXCESSDAYSLEEP INSOMNIA SLEEPDUR eduyears 
do
awk '{ if (NR!=1) print $2,$3}' $i > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/IVEA/"$i"EA
done

for i in blonde brown red 
do
awk '{ if (NR!=1) print $1,$2}' $i > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/IVEA/"$i"EA
done

#4. EXTRACT THE EFFECT ALLELE CODING FROM THE PD DATASER
cd /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/
mkdir /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/IVinPDEA/
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation EXCESSDAYSLEEP INSOMNIA SLEEPDUR eduyears blonde brown red
do
echo "RsID REF" > header
awk '{ if (NR!=1) print $10,$2}' "$i"SNPsinNallsPDdataV > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/IVinPDEA/"$i"inPDEA
done

#5.Find which snps Effect Allere are different between the trait and PD datasets, for then flipping them following the same approach as before (R script)
cd /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/
echo "CHRPOS A1 A2 freq b se p N_cases N_controls SNP" | sed -e 's/ /\t/g' > header
mkdir /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/SNPsinPD4Flipping/
mkdir /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/SameEAIVPD/

for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation EXCESSDAYSLEEP INSOMNIA SLEEPDUR eduyears blonde brown red
do
diff IVEA/"$i"EA IVinPDEA/"$i"inPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp ../"$i"SNPsinNallsPDdataV > SameEAIVPD/"$i"PDHM
grep -f temp ../"$i"SNPsinNallsPDdataV > prueba
cat header prueba > SNPsinPD4Flipping/"$i"inPD4flip
done

mkdir /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/PDSNPsflipped

#*inPD4flip will be the input for the second step of the Rscript titled Harmonization.Rmd, which will returned the *fliped files

#6.After fliping the SNPs that needed to be fliped, join them with the ones that already had the same Effect allele coding as their corresponding exposure GWAS data.
d /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/PDSNPsflipped/
mkdir /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/SNPsinPDHarmony/

for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation EXCESSDAYSLEEP INSOMNIA SLEEPDUR eduyears blonde brown red
do
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' $i > temp
cat /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/SameEAIVPD/"$i"PDHM temp | sed -e 's/ /\t/g'| sort > /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/SNPsinPDHarmony/"$i"PDinHarmony
done

#/data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/SNPsinPDHarmony/ will contain the PD gwas data that will serve as the input for MR
#/data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/IVinHarmonyPositive/ will contain each of the traits gwas data that will serve as the input for MR
