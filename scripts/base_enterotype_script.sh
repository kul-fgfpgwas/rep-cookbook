#/bin/bash

#PBS -N Enterotype
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00

### What is my phenotype to test ?
VAR1=enterotype

### NORMALIZATION
NORM=ENTERO

### Run by chromosome
CHR=data_chr21.bgen.gen

### move to directory with genotype data
cd /yourcompletedirectorypath/

## make a new directory, base on phenotype, to place my data in:
mkdir /yourcompletedirectorypath/gwasresults/${VAR1}
mkdir /yourcompletedirectorypath/gwasresults/${VAR1}/${NORM}
mkdir /yourcompletedirectorypath/gwasresults/${VAR1}/${NORM}/log/

### load my snptest module
module add apps/snptest.2.5.4-BETA3 

time snptest_v2.5.4-beta3 -data ${CHR} /yourcompletedirectorypath/FGFP.sample -o /yourcompletedirectorypath/gwasresults/$VAR1/${NORM}/${CHR}.out -log /yourcompletedirectorypath/$VAR1/${NORM}/log/${CHR}.log -frequentist add -method newml -pheno ${VAR1} -baseline_phenotype 1
