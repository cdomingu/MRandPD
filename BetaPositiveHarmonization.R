########################################### BetaPositiveHarmonization.R =====
# 
# author: Carmen Domínguez
# email: cdomingu@lcg.unam.mx
# date: 11/22/2019
# version: 1.0
# Usage: Rscript --vanilla /data/neurogen/MRandPD/Scripts/Rscripts/BetaPositiveHarmonization.R --input /data/neurogen/MRandPD/Results/IndexSNPs_inExposuresGWASdata/exposureIVs_inExposureGWASdata --exposure exposure
# This script is for ensuring that all betas are positive on the exposure dataset. In order to do so, it first identifies the betas that are already positive, and selects the betas that are negative.
# On the ones that are negative, it flips the effect and the other allele, in recalculates the effect allele frequency by substracting it from 1, and changes the sign of the beta by multiplying it by -1. 
# Lastly, it rejoins the data that was just harmonized with the one that already had a positive beta.

###########################################
#library("optparse")
#options(stringsAsFactors=F)

#option_list = list(
#  make_option(c("-i", "--input"), type="character", default=NULL,
#              help="Input expression dataset file name", metavar="character"),
#  make_option(c("-e", "--exposurename"), type="character", default=NULL,
#                help="Exposure of interest's name", metavar="character")
#);

#opt_parser = OptionParser(option_list=option_list);
#opt = parse_args(opt_parser);
#input_exposure4flip_filename=opt$input
#exposure_name=opt$exposurename

args <- commandArgs()
input_exposure_filename=args[8]
exposure_name=args[10]

IVcolnames<-c("RsID","chr","pos","EA","OA","EAF","beta","se","pvalue")
exposureIVsdata<-read.table(input_exposure_filename,skip=1, col.names=IVcolnames)
Betapos <- exposureIVsdata[exposureIVsdata$beta > 0, IVcolnames] 
Betaneg <- exposureIVsdata[exposureIVsdata$beta < 0, IVcolnames]
Betaneg[,4:5]<- Betaneg[,5:4]
Betaneg[,6]<- 1-Betaneg[,6]
Betaneg[,7]<- Betaneg[,7]*-1
new_exposureIVsdata<-rbind(Betapos,Betaneg)
write.table(new_exposureIVsdata,file =paste("/data/neurogen/MRandPD/Results/BetaPositiveHarmonized/",exposure_name,"_betapositive",sep = "") , sep="\t", quote = FALSE,row.names = F)
