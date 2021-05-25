#!/bin/Rscript

#This is a modified script from the original: https://github.com/explodecomputer/mr_frailty/blob/master/inst/analysis/bmi_pd/scripts/aggregate_model1.R
#from: Noyce, A. J. et al. Estimating the causal influence of body mass index on risk of Parkinson disease: A Mendelian randomisation study. PLoS Med. 14, 1–19 (2017).

library(dplyr)

splits <- 100
outdir <- "../FrailtyAnalysis/analysis/dpw_pd/dpwscratch"
savefile <- "../FrailtyAnalysis/analysis/dpw_pd/results/alldpwmodel.RData"

l <- list()
for(i in 1:splits)
{
	cat(i, "\n")
	filename <- paste0(outdir, "/modelresults", i, ".RData")
	if(file.exists(filename))
	{
		load(filename)
		l[[i]] <- res
	} else {
		message("Missing ", filename)
	}
}

res <- bind_rows(l)
save(res, file = savefile)
