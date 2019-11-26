########################################### R2andFstatistic_calculation.R =====
# 
# author: Carmen Domínguez
# email: cdomingu@lcg.unam.mx
# date: 11/26/2019
# version: 1.0
# Usage: Rscript --vanilla /data/neurogen/MRandPD/Scripts/Rscripts/R2andFstatistic_calculation.R --input /data/neurogen/MRandPD/Results/R2andFstatCalculation/exposureIVs_inExposureGWASdata --exposure exposure
# This script is for calculating the percetage of variance in the exposure explained by the Instrumental Variables (Rsquare)
# As well as calculatin the F statistic, a measure of the strenght of the association between the IVs and the exposure
###########################################
args <- commandArgs()
input_filename=args[8]
exposure_name=args[10]

#To get the Rsquare I followed the formula that I found on: https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0120758.s001 which is basically the same used by Gharahkhani,2019 (BMI&cancer):
#2(β^2)MAF(1 − MAF)/ 2(β^2)MAF(1 − MAF) + (se(β))^2*2NMAF(1 − MAF)

r2function <- function(x) {
  Rsq <- (2*(x[2]^2)*x[1]*(1-x[1]))/((2*(x[2]^2)*x[1]*(1-x[1]))+((x[3]^2)*2*x[4]*x[1]*(1-x[1])))
  result<- paste(result,Rsq)
  return(result)
}

#To calculate the F-statistic I followed the formula on Noyce,2017 (MBI&Parkinson). Rsquare(n-1-k)/(1-Rsquare)k
#where n--> sample size and k--> number of variants
#(According to Gharahkhani,2019 (BMI&cancer), "the genetic variants used as the instrumental variable must be strongly associated with the risk factor of interest (F-statistics > 10)".)

exposureforR2calc <-read.table(input_filename)
result<- vector()
RsqIV<- apply(exposureforR2calc[,c(6,7,8,10)],1,r2function)
Rsquare<- sum(as.numeric(RsqIV))  
PercentaceVarExplained <- Rsquare*100
Fstatistic <- Rsquare*(max(exposureforR2calc[,10])-1-nrow(exposureforR2calc))/((1-Rsquare)*nrow(exposureforR2calc))
cat("The percentage of variance explained is: ",PercentaceVarExplained)
cat("\n")
cat("The F statistic for this study is : ", Fstatistic)
cat("\n")
