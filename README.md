# rep-cookbook
Replication Cookbook for FGFP GWAS

## (i) Microbiome processing: Sample inference from amplicon data
Perform sample inference from high-throughput amplicon data using the SixteenS pipeline and [dada2](https://github.com/benjjneb/dada2).

### SixteenS pipeline performs the following tasks:
1. Merging reads: Forward and reverse reads are merged with the FLASH software (v1.2).
2. Reads QC: quality filtering was subsequently performed with the fastx toolkit (http://hannonlab.cshl.edu/fastx_toolkit/), excluding those sequences with >5% nucleotides of quality score <30.
3. Chimera removal: chimeras in sequences were removed using UCHIME (v6.0) with gold.fa as reference. 
4. Subsampling: randomly selecting 10,000 reads for each sample.
5. Taxonomical classification and compositional matrices for each taxonomical level were carried out using RDP classifier with RDP14, where classifications with low confidence at the genus level (<0.8) were organized in an arbitrary taxon of "unclassified_family".

**Run the analysis (in bash)**
```
bash SixteenS.sh
```

### dada2 pipeline performs the following steps:
1. Quality control: Perform filtering and trimming of reads.
2. Perform dereplication: Group amplicon reads with same sequence into unique sequences.
3. Sample inference: Using the core sequence-variant inference algorithm to the dereplicated data.
4. Remove chimeras: Identify sequences that are exact bimeras.
5. Taxonomical classification.

**Run the  analysis (in R):**
```
Seq_dada2.R
```

## (ii) Genotype QC and imputation
### Description
1. FGFP samples were Samples processed on two different arrays, Human Core Exome v1.0: 576 samples and Human Core Exome v1.1: 2112 samples, resulting in 545, 535 markers in merged raw data set. 
2. Allele calling done by GenomeStudio
3. SNP QC: cross check with sex information, remove Unmapped (chr0) Variants, update to RS ids and remove duplicates; remove variants with >5% missing, remove monomorphic alleles, remove HWE p<1e-5 and remove AT/CG ambiguous sites. 
4. Sample QC: remove sample with >5% SNPs missing, remove samples with heterozygosity <3SD of mean value, remove samples with relatedness >0.025. 
5. Merge with 1KG samples and remove those with outstanding population structure. 
6. Imputation with HRCv1.1 server and EUR population as reference. 
7. Post imputation QC: remove SNPs that are monomorphic, INFO < 0.5 & MAF < 1% .

### Method:
*needs to be added*

## (iii) Association tests
1. Select overlapping CMMs between FGFP and Kiel cohorts (CMM list “TaxaNamesAll.txt”).
   - There are a total of 72 taxa and 3 alpha diversity measures
2. Transformations and association analysis. We performed:
   - Linear regression in SNPTEST.
     - Log2 transformed data
     - Rank normal transformed data
    - Generalized linear regression on binary traits in SNPTEST.
      - This is performed on all potential zero-inflated taxa in a hurdle step analysis
        a. Hurdle taxa = those taxa with 5% of more individuals with zero counts
        b. Zero samples set to “0”, non-zero samples set to “1” to for a test of presence vs absence
    - Lude Franke meta analysis pipeline.
   The Rscript “fgfp_gwas_cookbook_help.R” is provided to assist in transformation of the taxa abundance and alpha-diversity data. Steps:
      1. Identify the taxa with >=5% of samples with zero-values. 
         a. While the above step would be the typical procedure we would like for you define the same taxa that we did as the hurdle taxa. These taxa are listed in the file “TaxaNamesHurdleOnly.txt”.
         b. Define these as hurdle taxa
         c. create a binary variable for these taxa
            i. 0 = 0
            ii. !0 =1
         d. create zero-truncated abundance data for these taxa
            i. all zero values are turned to NA and all non-zero values are treated as normal.
      2. All remaining abundance data (both non-hurdle and hurdle-truncated abundance data) are rank normal transformed and log2 transformed
         a. Rank normal transformation: completed with the rntransfrom() function from the GenABEL package
         b. Log2 transformation: all 0 values are turned into NA, effectively performing a truncation of zero values just as done in the hurdle step.
      3. Residualize RNT and LOG2 data
         a. Fit a linear model in R, or your favorite programing language, setting the abundance values as the response and the following variables as explanatory variables. 
            i. Top 10 PCs
            ii. Sex, age
            iii. Study specific batch variables, such as plate, processing date, etc..
         b. Extract the residuals of the fitted model
         c. Take care to return NAs for those samples that were already NAs. This is done to maintain order and structure in the data file for the association analysis

3. Perform the association in SNPTEST:
   a. Construct the .sample file for SNPTEST that includes all of the 
      i. Sample ids
      ii. Covariates
      iii. LOG2 residual data
      iv. RNT residual data
      v. Hurdle Binary data
   b. Extract Variants to test
      i. The ~25K variant positions to test are in the file “Sites_2_Test.txt”
      ii. This file contains HRCv1.1 SNPs and 1000G indels, as such you may not identify all of the variants in your data. This is OK! Proceed with all of those that you can identify.
          1. We find 23,098 significant variants in the HRCv1.1 data.
      iii. Extract and prepare genotype data with qctools
           1. Below is some code that could be used to do that

```
for i in ~/YOURDIRECTORYPATH/*.bgen; do echo ${i}; qctool -g ${i} -incl-rsids SNPids2test.txt -ofiletype gen -og sigsites_chr${i# YOURDIRECTORYPATH /data_chr}.gen; done 

wc -l sigsites_chr*.gen  # 23,098 sites identified

## Concatenate all the sigsites together
for i in sigsites_chr*.gen; do echo ${i}; cat ${i} >> all_sigsites2text.gen; done
wc -l all_sigsites2text.gen # 23,098 sites 

## remove chromosome data
rm sigsites_chr*.gen
```
   c. Perform association analysis for each taxa transformation
      i. EXAMPLES

```
RUN SNPTEST for Hurdle Binary Analysis

VAR1=C_Alphaproteobacteria_HB 
snptest_v2.5 -data all_sigsites2text.gen TransResCovar_MYSTUDY.sample -o $VAR1.out -log $VAR1.log -frequentist 1 -method score -cov_names PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 batchvar1 batchvar2 sex age -pheno $VAR1 -use_raw_phenotypes


RUN SNPTEST for linear RNT and LOG analysis

VAR1=G_Blautia_RNT_Residuals
snptest_v2.5 -data all_sigsites2text.gen TransResCovar_MYSTUDY.sample -o $VAR1.out -log $VAR1.log -frequentist 1 -method score -pheno $VAR1 -use_raw_phenotypes
```

