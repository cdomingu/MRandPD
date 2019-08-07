#This Script is for conducting Harmonization to ensure that both datasets that contain the indexSNPs for each exposure (IV-exposure, IV-outcome),
#are identically coded regarding the effect allele. If not, we will flip the one from the PD GWAS Dataset so that it matches its exposure counterpart

#1.Convert all the alleles to Upper case, as in the PD GWAS dataset they are in lowercase
cd /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/
sed 's/\([a-g]\)/\U\1/g' cigperdaySNPinPDdataV > temp
sed 's/\([t]\)/\U\1/g' temp | sort -k1 > cigperdaySNPinPDdataV
sed 's/\([a-g]\)/\U\1/g' ageofinitsmkSNPinPDdataV > temp
sed 's/\([t]\)/\U\1/g' temp | sort -k1 > ageofinitsmkSNPinPDdataV
sed 's/\([a-g]\)/\U\1/g' fulldrinksperweekSNPinPDdataV > temp
sed 's/\([t]\)/\U\1/g' temp | sort -k1 > fulldrinksperweekSNPinPDdataV
sed 's/\([a-g]\)/\U\1/g' fulleduyearsSNPinPDdataV > temp
sed 's/\([t]\)/\U\1/g' temp | sort -k1 > fulleduyearsSNPinPDdataV
sed 's/\([a-g]\)/\U\1/g' fulleversmkSNPinPDdataV > temp
sed 's/\([t]\)/\U\1/g' temp | sort -k1 > fulleversmkSNPinPDdataV
sed 's/\([a-g]\)/\U\1/g' fullexcessdaysleepSNPinPDV > temp
sed 's/\([t]\)/\U\1/g' temp | sort -k1 > fullexcessdaysleepSNPinPDV
sed 's/\([a-g]\)/\U\1/g' insomniaSNPinPDdataV > temp
sed 's/\([t]\)/\U\1/g' temp | sort -k1 > insomniaSNPinPDdataV
sed 's/\([a-g]\)/\U\1/g' sleepdurSNPinPDdataV > temp
sed 's/\([t]\)/\U\1/g' temp | sort -k1 > sleepdurSNPinPDdataV
sed 's/\([a-g]\)/\U\1/g' smokingsatSNPinPDdataV > temp
sed 's/\([t]\)/\U\1/g' temp | sort -k1 > smokingsatSNPinPDdataV

#2. Extract only the arkerID, the Reference allele and the alternative allele from the SNPs present in the PD GWAS data for each exposure
echo "RsID REF ALT" > header
awk '{ if (NR!=1) print $1,$2,$3}' cigperdaySNPinPDdataV > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/cigarettespdayinPDEA
awk '{ if (NR!=1) print $1,$2,$3}' ageofinitsmkSNPinPDdataV > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/ageofinitsmokinginPDEA
awk '{ if (NR!=1) print $1,$2,$3}' fulldrinksperweekSNPinPDdataV > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/drinksperweekinPDEA
awk '{ if (NR!=1) print $1,$2,$3}' fulleduyearsSNPinPDdataV > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/eduyearsinPDEA
awk '{ if (NR!=1) print $1,$2,$3}' fulleversmkSNPinPDdataV > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/eversmokinginPDEA
awk '{ if (NR!=1) print $1,$2,$3}' fullexcessdaysleepSNPinPDV > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/excessdaysleepinPDEA
awk '{ if (NR!=1) print $1,$2,$3}' insomniaSNPinPDdataV > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/insomniainPDEA
awk '{ if (NR!=1) print $1,$2,$3}' sleepdurSNPinPDdataV > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/sleepdurinPDEA
awk '{ if (NR!=1) print $1,$2,$3}' smokingsatSNPinPDdataV > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/smokingstatusinPDEA

#3. Extract only the arkerID, the Reference allele and the alternative allele from the SNPs present in the PD GWAS data for each exposure
cd /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/
echo "RsID REF ALT" > header
awk '{ if (NR!=1) print $2,$3,$4}' cigarettespday > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/cigarettespdayEA
awk '{ if (NR!=1) print $2,$3,$4}' ageofinitsmoking > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/ageofinitsmokingEA
awk '{ if (NR!=1) print $2,$3,$4}' drinksperweek > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/drinksperweekEA
awk '{ if (NR!=1) print $2,$3,$4}' eversmoking > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/eversmokingEA
awk '{ if (NR!=1) print $2,$3,$4}' smokingstatus > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/smokingstatusEA
awk '{ if (NR!=1) print $2,$3,$4}' eduyears > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/eduyearsEA
awk '{ if (NR!=1) print $2,$3,$4}' excessdaysleep > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/excessdaysleepEA
awk '{ if (NR!=1) print $2,$3,$4}' insomnia > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/insomniaEA
awk '{ if (NR!=1) print $2,$3,$4}' sleepdur > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/sleepdurEA

#4.Find which snps Effect Allere are different between the trait and PD datasets, for then flipping them following the same approach as before (R script)
cd /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/
diff ageofinitsmokingEA ageofinitsmokinginPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/ageofinitsmkSNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/AllHarmonized/ageofinitPDHM
grep -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/ageofinitsmkSNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPD4Flipping/ageofinitPD4flip
diff cigarettespdayEA cigarettespdayinPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/cigperdaySNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/AllHarmonized/cigarettespdayinPDHM
grep -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/cigperdaySNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPD4Flipping/cigarettespdayinPD4flip
diff drinksperweekEA drinksperweekinPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/fulldrinksperweekSNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/AllHarmonized/drinksperweekinPDHM
grep -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/fulldrinksperweekSNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPD4Flipping/drinksperweekinPD4flip
diff eversmokingEA eversmokinginPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/fulleversmkSNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/AllHarmonized/eversmokinginPDHM
grep -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/fulleversmkSNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPD4Flipping/eversmokinginPD4flip
diff smokingstatusEA smokingstatusinPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/smokingsatSNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/AllHarmonized/smokingstatusinPDHM
grep -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/smokingsatSNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPD4Flipping/smokingstatusinPD4flip
diff eduyearsEA eduyearsinPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/fulleduyearsSNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/AllHarmonized/eduyearsinPDHM
grep -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/fulleduyearsSNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPD4Flipping/eduyearsinPD4flip
diff insomniaEA insomniainPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/insomniaSNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPD4Flipping/insomniainPD4flip
diff sleepdurEA sleepdurinPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/sleepdurSNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/AllHarmonized/sleepdurinPDHM
grep -f temp /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/sleepdurSNPinPDdataV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPD4Flipping/sleepdurinPD4flip
grep -f excessdaysleepinPDEA /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/SNPsinPDdata/fullexcessdaysleepSNPinPDV > /data/neurogen/MRandPD/Results/HarmonizationVerbose/AllHarmonized/excessdaysleepPDinHarmony
#*inPD4flip will be the input for the Rscript titled HarmonSNPinPDdata.Rmd, which will returned the *fliped files

#5.After fliping the SNPs that needed to be fliped, join them with the ones that already had the same Effect allele coding as their corresponding exposure GWAS data.
cd /data/neurogen/MRandPD/Results/HarmonizationVerbose/data4Harmonizing/AllHarmonized/
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9}' ageofinitPDfliped | sed -e 's/"//g' | sed -e 's/ /\t/g' > temp
cat ageofinitPDHM temp | sort > ageofinitPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9}' cigarettespdayPDfliped | sed -e 's/"//g' | sed -e 's/ /\t/g' > temp
cat cigarettespdayinPDHM temp | sort > cigarettespdayPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9}' drinksperweekPDfliped | sed -e 's/"//g' | sed -e 's/ /\t/g' > temp
cat drinksperweekinPDHM temp | sort > drinksperweekPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9}' eduyearsPDfliped | sed -e 's/"//g' | sed -e 's/ /\t/g' > temp
cat eduyearsinPDHM temp | sort > eduyearsPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9}' eversmokingPDfliped | sed -e 's/"//g' | sed -e 's/ /\t/g' > temp
cat eversmokinginPDHM temp | sort > eversmokingPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9}' sleepdurPDfliped | sed -e 's/"//g' | sed -e 's/ /\t/g' > temp
cat sleepdurinPDHM temp | sort > sleepdurPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9}' smokingstatusPDfliped | sed -e 's/"//g' | sed -e 's/ /\t/g' > temp
cat smokingstatusinPDHM temp | sort > smokingstatusPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9}' insomniaPDfliped | sed -e 's/"//g' | sed -e 's/ /\t/g' > insomniaPDinHarmony
