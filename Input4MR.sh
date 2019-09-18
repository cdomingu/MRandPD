#This script is for creating the INPUT files for each of the MR methods: IVW, MRPRESSO and GSMR
cd /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/IVinHarmonyPositive/
mkdir /data/neurogen/MRandPD/Results/Input4IVW/
mkdir /data/neurogen/MRandPD/Results/Input4MRPRESSO/
mkdir /data/neurogen/MRandPD/Results/Input4GSMR/

echo "bx bxse by byse snps ea oa eaf" > headerIVW
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation drinksperweek EXCESSDAYSLEEP INSOMNIA SLEEPDUR eduyears 
do
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7}' $i | sort > tempIV
awk '{ if (NR!=1) print $5,$6,$10}' /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/SNPsinPDHarmony/"$i"PDinHarmony | sort -k3 > tempPD
join -1 1 -2 3 tempIV tempPD > temp
awk '{ print $5,$6,$7,$8,$1,$2,$3,$4}' temp > PDIV
cat headerIVW PDIV | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/Input4IVW/"$i"IVWinput
done 

echo "bx bxse bxpval by byse bypval" > headerMRPRESSO
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation drinksperweek EXCESSDAYSLEEP INSOMNIA SLEEPDUR eduyears
do
awk '{ if (NR!=1) print $2,$6,$7,$8}' $i | sort > tempIV
awk '{ if (NR!=1) print $5,$6,$7,$10}' /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/SNPsinPDHarmony/"$i"PDinHarmony | sort -k4 > tempPD
join -1 1 -2 4 tempIV tempPD > temp
awk '{print $2,$3,$4,$5,$6,$7}' temp > PDIV
cat headerMRPRESSO 	PDIV | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/Input4MRPRESSO/"$i"MRPRESSOinput
done

echo "SNP a1 a2 a1freq bx bxse bxpval bxn by byse bypval byn" > headerGSMR
for i in CigarettesPerDay EverSmoking SmokingStatus AgeOfInitiation drinksperweek EXCESSDAYSLEEP INSOMNIA SLEEPDUR eduyears
do
awk '{ if (NR!=1) print $2,$3,$4,$5,$6,$7,$8,$9}' $i | sort > tempIV
awk '{ if (NR!=1) print $5,$6,$7,$8,$10}' /data/neurogen/MRandPD/Results/HarmonizationNallsPDdata/data4Harmony/SNPsinPDHarmony/"$i"PDinHarmony | sort -k5> tempPD
join -1 1 -2 5 tempIV tempPD > IVPD
cat headerGSMR IVPD | sed -e 's/ /\t/g' > /data/neurogen/MRandPD/Results/Input4GSMR/"$i"GSMRinput
done
