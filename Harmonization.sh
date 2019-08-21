#This Script is for conducting Harmonization to ensure that both datasets that contain the indexSNPs for each exposure (IV-exposure, IV-outcome),
#are identically coded regarding the effect allele. If not, we will flip the one from the PD GWAS Dataset so that it matches its exposure counterpart


#2. Extract the information from the SNPs present in the PD GWAS data for each exposure
cd /data/neurogen/MRandPD/Results/R2andFstatisticVerbose/
echo "RsID REF ALT AF BETA SE PVAL" > header
awk '{ if (NR!=1) print $1,$4,$5,$6,$7,$8,$9}' fullexcessdaysleepforR2calcVNalls > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/InstruVar/excessdaysleepIV
awk '{ if (NR!=1) print $1,$4,$5,$6,$7,$8,$9}' sleepdurforR2calcVNalls > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/InstruVar/sleepdurIV
awk '{ if (NR!=1) print $1,$4,$5,$6,$7,$8,$9}' insomniaforR2calcVNalls > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/InstruVar/insomniaIV
awk '{ if (NR!=1) print $1,$4,$5,$6,$7,$8,$9}' eduyearsforR2calcVNalls > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/InstruVar/eduyearsIV

cd /data/neurogen/MRandPD/Results/R2andFstatisticprePD/
echo "RsID REF ALT AF BETA SE PVAL" > header
awk '{ if (NR!=1) print $3,$4,$5,$6,$9,$10,$8}' ageofinitiationforR2calc > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/InstruVar/agepfinitiationIV
awk '{ if (NR!=1) print $3,$4,$5,$6,$9,$10,$8}' cigarettespdayforR2calc > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/InstruVar/cigarettespdayIV
awk '{ if (NR!=1) print $3,$4,$5,$6,$9,$10,$8}' drinksperweekforR2calc > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/InstruVar/drinksperweekIV
awk '{ if (NR!=1) print $3,$4,$5,$6,$9,$10,$8}' eversmokingforR2calc > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/InstruVar/eversmokingIV
awk '{ if (NR!=1) print $3,$4,$5,$6,$9,$10,$8}' smokingstatusforR2calc > temp
cat header temp | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/InstruVar/smokingstatusIV

#3. AFTER HARMONIZING THE IV TO ALL POSITIVE EFFECT SIZE (Harmonizarion.rmd) AND EXTRACT THE EFFECT ALLLELE CODING
cd /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/IVinHarmonyPositive/
echo "RsID REF ALT" > header
awk '{ if (NR!=1) print $2,$3,$4}' cigarettespday > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/IVEA/cigarettespdayEA
awk '{ if (NR!=1) print $2,$3,$4}' ageofinitsmoking > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/IVEA/ageofinitsmokingEA
awk '{ if (NR!=1) print $2,$3,$4}' drinksperweek > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/IVEA/drinksperweekEA
awk '{ if (NR!=1) print $2,$3,$4}' eversmoking > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/IVEA/eversmokingEA
awk '{ if (NR!=1) print $2,$3,$4}' smokingstatus > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/IVEA/smokingstatusEA
awk '{ if (NR!=1) print $2,$3,$4}' eduyears > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/IVEA/eduyearsEA
awk '{ if (NR!=1) print $2,$3,$4}' excessdaysleep > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/IVEA/excessdaysleepEA
awk '{ if (NR!=1) print $2,$3,$4}' insomnia > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/IVEA/insomniaEA
awk '{ if (NR!=1) print $2,$3,$4}' sleepdur > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/IVEA/sleepdurEA

#5. EXTRACT THE EFFECT ALLELE CODING FROM THE PD DATASER
cd /data/neurogen/MRandPD/Results/HarmonizationNalls2019/
echo "RsID REF ALT" > header
awk '{ if (NR!=1) print $10,$2,$3}' cigarettespdaySNPinNallsPDdataV > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > data4Harmony/IVinPDEA/cigarettespdayinPDEA
awk '{ if (NR!=1) print $10,$2,$3}' ageofinitSNPinNallsPDdataV > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > data4Harmony/IVinPDEA/ageofinitsmokinginPDEA
awk '{ if (NR!=1) print $10,$2,$3}' drinkspweekSNPinNallsPDdataV > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > data4Harmony/IVinPDEA/drinksperweekinPDEA
awk '{ if (NR!=1) print $10,$2,$3}' eversmokingSNPinNallsPDdataV > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > data4Harmony/IVinPDEA/eversmokinginPDEA
awk '{ if (NR!=1) print $10,$2,$3}' smokingstatSNPinNallsPDdataV > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > data4Harmony/IVinPDEA/smokingstatusinPDEA
awk '{ if (NR!=1) print $10,$2,$3}' fullexcessdaysleepSNPinNallsPDdataV > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > data4Harmony/IVinPDEA/excessdaysleepinPDEA
awk '{ if (NR!=1) print $10,$2,$3}' insomniaSNPinNallsPDdataV > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > data4Harmony/IVinPDEA/insomniainPDEA
awk '{ if (NR!=1) print $10,$2,$3}' sleepdurSNPinNallsPDdataV > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > data4Harmony/IVinPDEA/sleepdurinPDEA
awk '{ if (NR!=1) print $10,$2,$3}' fulleduyearsSNPinNallsPDdataV > temp
cat header temp | sed -e 's/ /\t/g' | sort -k1 > data4Harmony/IVinPDEA/eduyearsinPDEA

#5.Find which snps Effect Allere are different between the trait and PD datasets, for then flipping them following the same approach as before (R script)
echo "CHRPOS A1 A2 freq b se p N_cases N_controls SNP" | sed -e 's/ /\t/g' > header
cd /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/
diff IVEA/ageofinitsmokingEA IVinPDEA/ageofinitsmokinginPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp ../ageofinitSNPinNallsPDdataV > SameEAIVPD/ageofinitPDHM
grep -f temp ../ageofinitSNPinNallsPDdataV > prueba
cat header prueba > SNPsinPD4Flipping/ageofinitinPD4flip
diff IVEA/cigarettespdayEA IVinPDEA/cigarettespdayinPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp ../cigarettespdaySNPinNallsPDdataV > SameEAIVPD/cigarettespdayinPDHM
grep -f temp ../cigarettespdaySNPinNallsPDdataV > prueba
cat header prueba > SNPsinPD4Flipping/cigarettespdayinPD4flip
diff IVEA/drinksperweekEA IVinPDEA/drinksperweekinPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp ../drinkspweekSNPinNallsPDdataV > SameEAIVPD/drinksperweekinPDHM
grep -f temp ../drinkspweekSNPinNallsPDdataV > prueba
cat header prueba > SNPsinPD4Flipping/drinksperweekinPD4flip
diff IVEA/eversmokingEA IVinPDEA/eversmokinginPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp ../eversmokingSNPinNallsPDdataV > SameEAIVPD/eversmokinginPDHM
grep -f temp ../eversmokingSNPinNallsPDdataV > prueba
cat header prueba > SNPsinPD4Flipping/eversmokinginPD4flip
diff IVEA/smokingstatusEA IVinPDEA/smokingstatusinPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp ../smokingstatSNPinNallsPDdataV > SameEAIVPD/smokingstatusinPDHM
grep -f temp ../smokingstatSNPinNallsPDdataV > prueba
cat header prueba > SNPsinPD4Flipping/smokingstatusinPD4flip
diff IVEA/eduyearsEA IVinPDEA/eduyearsinPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp ../fulleduyearsSNPinNallsPDdataV > SameEAIVPD/eduyearsinPDHM
grep -f temp ../fulleduyearsSNPinNallsPDdataV > prueba
cat header prueba > SNPsinPD4Flipping/eduyearsinPD4flip
diff IVEA/insomniaEA IVinPDEA/insomniainPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp ../insomniaSNPinNallsPDdataV > SameEAIVPD/insomniainPDHM
grep -f temp ../insomniaSNPinNallsPDdataV > prueba
cat header prueba > SNPsinPD4Flipping/insomniainPD4flip
diff IVEA/sleepdurEA IVinPDEA/sleepdurinPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp ../sleepdurSNPinNallsPDdataV > SameEAIVPD/sleepdurinPDHM
grep -f temp ../sleepdurSNPinNallsPDdataV > prueba
cat header prueba > SNPsinPD4Flipping/sleepdurinPD4flip
diff IVEA/excessdaysleepEA IVinPDEA/excessdaysleepinPDEA | sed -e 's/ /\t/g' | cut -f2 | grep 'rs' | sort | uniq > temp
grep -v -f temp ../fullexcessdaysleepSNPinNallsPDdataV > SameEAIVPD/excessdaysleepinPDHM
grep -f temp ../fullexcessdaysleepSNPinNallsPDdataV > prueba 
cat header prueba > SNPsinPD4Flipping/excessdaysleepinPD4flip

#*inPD4flip will be the input for the Rscript titled HarmonSNPinPDdata.Rmd, which will returned the *fliped files

#5.After fliping the SNPs that needed to be fliped, join them with the ones that already had the same Effect allele coding as their corresponding exposure GWAS data.
cd /data/neurogen/MRandPD/Results/HarmonizationNalls2019/data4Harmony/PDSNPsflipped
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' ageofinitPDfliped > temp
cat ../SameEAIVPD/ageofinitPDHM temp | sed -e 's/ /\t/g'| sort > ../SNPsinPDHarmony/ageofinitPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' cigarettespdayPDfliped > temp
cat ../SameEAIVPD/cigarettespdayinPDHM temp | sed -e 's/ /\t/g'| sort > ../SNPsinPDHarmony/cigarettespdayPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' drinksperweekPDfliped > temp
cat ../SameEAIVPD/drinksperweekinPDHM temp | sed -e 's/ /\t/g'| sort > ../SNPsinPDHarmony/drinksperweekPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' eduyearsPDfliped > temp
cat ../SameEAIVPD/eduyearsinPDHM temp | sed -e 's/ /\t/g'| sort > ../SNPsinPDHarmony/eduyearsPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' eversmokingPDfliped > temp
cat ../SameEAIVPD/eversmokinginPDHM temp | sed -e 's/ /\t/g'| sort > ../SNPsinPDHarmony/eversmokingPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' sleepdurPDfliped > temp
cat ../SameEAIVPD/sleepdurinPDHM temp | sed -e 's/ /\t/g'| sort > ../SNPsinPDHarmony/sleepdurPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' smokingstatusPDfliped > temp
cat ../SameEAIVPD/smokingstatusinPDHM temp | sed -e 's/ /\t/g'| sort > ../SNPsinPDHarmony/smokingstatusPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' insomniaPDfliped > temp
cat ../SameEAIVPD/insomniainPDHM temp | sed -e 's/ /\t/g'| sort > ../SNPsinPDHarmony/insomniaPDinHarmony
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' excessdaysleepPDfliped > temp
cat ../SameEAIVPD/excessdaysleepinPDHM temp | sed -e 's/ /\t/g'| sort > ../SNPsinPDHarmony/excessdaysleepPDinHarmony
