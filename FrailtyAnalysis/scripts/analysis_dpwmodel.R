#!/bin/Rscript

#This is a modified script from the original: https://github.com/explodecomputer/mr_frailty/blob/master/inst/analysis/bmi_pd/scripts/analysis.R
#from: Noyce, A. J. et al. Estimating the causal influence of body mass index on risk of Parkinson disease: A Mendelian randomisation study. PLoS Med. 14, 1–19 (2017).

## ---- setup ----

library(ggplot2)
library(dplyr)
library(systemfit)
library(Runuran)

load("/export/storage/users/cdomingu/FrailtyAnalysis/analysis/dpw_pd/results/alldpwmodel.RData")
res$model <- "model1"
res$modname <- "DPW ~ SNP"


res_plot <- subset(res, test %in% c("2sls", "obs", "MR Egger", "Inverse variance weighted"))

res_plot$test <- as.factor(res_plot$test)
levels(res_plot$test) <- c("2SLS", "IVW", "MR Egger", "Obs assoc")


## ---- empirical_results ----

dat_ci <- dplyr::group_by(res_plot, test, model) %>% dplyr::summarise(b=mean(beta), lower_ci=quantile(beta, 0.05), upper_ci=quantile(beta, 0.95))

emp_dat <- data.frame(
	b = c(0.79, 0.73, 0.78),
	lower_ci = c(0.65, 0.56, 0.67),
	upper_ci = c(0.96, 0.95, 0.92),
	test = c("IVW", "MR Egger", "Obs assoc"),
	model = "Empirical"
)

dat_ci <- bind_rows(subset(dat_ci, test != "2SLS"), emp_dat)
dat_ci$model <- as.factor(dat_ci$model)
levels(dat_ci$model) <- c("True effect estimate,\nobtained from data", "Survival bias effect where \nDrinks per week \nhas no causal effect on PD,\nbased on 1000 simulations")
dat_ci
write.table(dat_ci,file="../FrailtyAnalysis/analysis/dpw_pd/images/dpw_pd2019reuslts.txt", sep =" ", quote= FALSE)
ggplot(dat_ci, aes(y=b, x=model,color=dat_ci$model)) +
geom_point() +
geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0) +
facet_grid(test ~ .) +
coord_flip() +
scale_color_manual(values=c("steelblue","firebrick3")) +
theme_bw() +
labs(x="", y="Odds ratio and 95% confidence intervals") +
theme(legend.position = "none")
ggsave("../FrailtyAnalysis/analysis/dpw_pd/images/dpw_pd2019.pdf")
