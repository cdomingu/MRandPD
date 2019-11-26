########################################### SameEAHarmonization.R =====
# 
# author: Carmen Domínguez
# email: cdomingu@lcg.unam.mx
# date: 11/22/2019
# version: 1.0
# Usage: Rscript --vanilla /data/neurogen/MRandPD/Scripts/Rscripts/SameEAHarmonization.R --input /data/neurogen/MRandPD/Results/EAHarmonizationNallsPDdata/SNPsinPD4Flipping/exposureIVs_inParkinsonsGWASdata_EA4flip --exposure exposure
# This script is for flipping the effect allele and its estimates from the variants on the Parkinson's data that have been shown to have a different allele coding from the exposure's dataset
# To do so, it flips the effect and the other allele, in recalculates the effect allele frequency by substracting it from 1, and changes the sign of the beta by multiplying it by -1. 

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
input_exposure4flip_filename=args[8]
exposure_name=args[10]

IVcolnames<-c("RsID","chr","pos","EA","OA","EAF","beta","se","pvalue","N")
exposureIVs_PDdata<-read.table(input_exposure4flip_filename, col.names=IVcolnames)
exposureIVs_PDdata[,4:5]<- exposureIVs_PDdata[,5:4]
exposureIVs_PDdata[,6]<- 1-exposureIVs_PDdata[,6]
exposureIVs_PDdata[,7]<- exposureIVs_PDdata[,7]*-1
write.table(exposureIVs_PDdata,file =paste("/data/neurogen/MRandPD/Results/SameEAHarmonized/",exposure_name,"_EAflipped",sep = ""), sep="\t", quote = FALSE,row.names = F, col.names = F)
