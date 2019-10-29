#!/bin/bash

#PBS -N snptest
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00

### What is my phenotype to test ?
VAR1=G_Akkermansia

### What is my normalization type ?
NORM=HB

### Run by chromosome
CHR=allreplicationmarkers.gen

### move to directory with genotype data
cd /yourcompletedirectorypath/

## make a new directory, base on phenotype, to place my data in:
mkdir /yourcompletedirectorypath/gwasresults/${VAR1}
mkdir /yourcompletedirectorypath/gwasresults/${VAR1}/${NORM}
mkdir /yourcompletedirectorypath/gwasresults/${VAR1}/${NORM}/log/

### load my snptest module
module add apps/snptest.2.5.0 

### Running a Chromosome at a time
time snptest_v2.5 -data ${CHR} /yourcompletedirectorypath/FGFP.sample -o /yourcompletedirectorypath/gwasresults/$VAR1/${NORM}/${CHR}.out -log /yourcompletedirectorypath/gwasresults/$VAR1/${NORM}/log/${CHR}.log -frequentist 1 -method score -cov_names PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 exttype extractyear aliquotyear aliquotby plate snpsex age -pheno ${VAR1}_${NORM} -use_raw_phenotypes
