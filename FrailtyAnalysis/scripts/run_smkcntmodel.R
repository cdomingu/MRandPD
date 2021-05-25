#!/bin/bash

#$ -N smkcnt_pd
#$ -o job_reports/smkcntmodel.out
#$ -e job_reports/smkcntmodel.err
#$ -l h_rt=4:00:00
#$ -t 1-100
#$ -S /bin/bash
#$ -pe smp 6

#This is a modified script from the original: https://github.com/explodecomputer/mr_frailty/blob/master/inst/analysis/bmi_pd/scripts/run_model1.sh
#from: Noyce, A. J. et al. Estimating the causal influence of body mass index on risk of Parkinson disease: A Mendelian randomisation study. PLoS Med. 14, 1–19 (2017).

set -e

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  SGE_TASK_ID=${1}
fi

i=${SGE_TASK_ID}
splits=10
outdir="../FrailtyAnalysis/analysis/smkcnt_pd/smkcntscratch"

Rscript --no-save --args ${i} ${splits} ${outdir} < ../FrailtyAnalysis/scripts/smkcntmodel.R
