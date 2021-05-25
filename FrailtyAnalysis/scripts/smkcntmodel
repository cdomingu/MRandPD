#!/bin/Rscript

#This is a modified script from the original: https://github.com/explodecomputer/mr_frailty/blob/master/inst/analysis/bmi_pd/scripts/model4.R
#from: Noyce, A. J. et al. Estimating the causal influence of body mass index on risk of Parkinson disease: A Mendelian randomisation study. PLoS Med. 14, 1–19 (2017).

#Model description
# exposure, outcome and mortality

# exposure has no effect on outcome

# exposure has an effect on mortality

# outcome occurs based on age

# mortality is a function of age and exposure



# exposure = snp(s)

# mortality = age + exposure

# outcome = age


# simulate population with snp(s)

# give them ages from a specified age distribution


# total_n
# alive_n

# snp_freq
# snp_beta
# snp_effect_allele



source("../FrailtyAnalysis/scripts/functions.R")
library(dplyr)
library(stringr)
library(systemfit)
library(TwoSampleMR)

main <- function()
{
	demog <- read.csv("../FrailtyAnalysis/data/pd2019_demographics.csv")
	age_summary <- get_age_summary(demog, "Cases", "Controls", "Case_age_mean", "Control_age_mean", "Case_age_sd", "Control_age_sd")
	smkcnt_snps <- read.table("../FrailtyAnalysis/data/smkcnt_Liu2019_clumped.txt", he=T)
	smkcnt_snps$b <- smkcnt_snps$b


	parameters <- expand.grid(sim = 1:1000)

	# Parallel
	arguments <- commandArgs(trailingOnly=TRUE)
        jid <- "all"
        outdir <- "../FrailtyAnalysis/analysis/smkcnt_pd/"
        if(length(arguments) > 0)
        {
                arguments <- unlist(str_split(arguments,"="))
                jid <- as.numeric(arguments[2])
                message("jid value is ", jid)
                splits <- as.numeric(arguments[4])
                message("splits value is ", splits)
                outdir <- arguments[5]
                message("outdir is ", outdir)
		stopifnot(all(!is.na(jid), !is.na(splits), !is.na(outdir)))

		first <- (jid - 1) * splits + 1
		last <- min(nrow(parameters), jid * splits)
		parameters <- parameters[first:last, , drop=FALSE]
	}

	# Set output file
	outfile <- paste(outdir, "/modelresults", jid, ".RData", sep="")
	message("Running ", jid, ": ", nrow(parameters), " rows")
	message("Saving in ", outfile)


	l1 <- list()
	l2 <- list()
	l3 <- list()
	l4 <- list()
	for (i in 1:nrow(parameters))
	{
		message(i)

		dat <- simulate_ages(age_summary$gn[3], age_summary$gm[3], age_summary$gs[3], max_age=100, min_age=40, sample_size_multiplier=4)
		dat$cc <- simulate_events(dat$age, NULL, pd_incidence)
		snps <- simulate_snps(nrow(dat), smkcnt_snps$a1freq)
		dat$smkcnt <- rbinom(nrow(dat), 1, 0.255)
		dat$alive <- simulate_events(dat$age, dat$smkcnt, smkcnt_survival)
		dat$grs <- snps %*% smkcnt_snps$b
		index <- sample_cases_controls(dat, age_summary, min_age=40, max_age=100)

		a <- summary(glm(cc ~ smkcnt, dat[index,], family="binomial"))
		b <- summary(glm(cc ~ grs, dat[index,], family="binomial"))
		c <- summary(systemfit(cc ~ smkcnt, "2SLS", inst = ~ grs, data = dat[index,]))

		l1[[i]] <- coefficients(a)[2,]
		l2[[i]] <- coefficients(b)[2,]
		l3[[i]] <- coefficients(c)[2,]

		ss <- get_summary_stats(dat, snps, index)
		mres <- do_mr(smkcnt_snps$b, ss$b, smkcnt_snps$se, ss$se)
		l4[[i]] <- data.frame(beta = mres$b, se = mres$se, tval = NA, pval = mres$pval, sim = parameters$sim[i], test = mres$method)
	}

	l1 <- as.data.frame(do.call(rbind, l1))
	names(l1) <- c("beta", "se", "tval", "pval")
	l1$sim <- parameters$sim
	l1$test <- "obs"

	l2 <- as.data.frame(do.call(rbind, l2))
	names(l2) <- c("beta", "se", "tval", "pval")
	l2$sim <- parameters$sim
	l2$test <- "grs"

	l3 <- as.data.frame(do.call(rbind, l3))
	names(l3) <- c("beta", "se", "tval", "pval")
	l3$sim <- parameters$sim
	l3$test <- "2sls"

	l4 <- bind_rows(l4)

	res <- rbind(l1, l2, l3, l4)

	save(res, file=outfile)

}




pd_incidence <- function(age, ...)
{
	pd_inc <- read.table("../FrailtyAnalysis/data/pd_incidence.txt", he=T)
	pd_inc$age <- (pd_inc$age_low + pd_inc$age_upp) / 2
	pd_inc$n <- pd_inc$personyears / ((pd_inc$age_upp + pd_inc$age_low) / 2)
	pd_inc$p <- pd_inc$cases / pd_inc$n
	mod <- lm(p ~ poly(age,4), weights=n, pd_inc)
	pred <- predict(mod, data.frame(age=age))
	pred[pred < 0] <- 0
	return(pred)
}


smkcnt_survival <- function(age, smkcnt)
{
	# hr of 1.62 for current smokers, based on Qiao et al., 2000
 	hr <- rep(1, length(smkcnt))
 	hr[smkcnt == 1] <- 1.62
 	survival <- (1 - gompertz_makeham_cdf(age)) ^ hr
 	return(survival)
}

main()
