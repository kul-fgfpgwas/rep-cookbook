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
#### 1 - Select overlapping CMMs between FGFP and Kiel cohorts.
  CMM list “TaxaNamesAll.txt”. There are a total of 72 taxa and 3 alpha diversity measures.

#### 2 - Transformations and association analysis.
   We performed:
   - Linear regression in SNPTEST.
     - Log2 transformed data
     - Rank normal transformed data
    - Generalized linear regression on binary traits in SNPTEST.
      - This is performed on all potential zero-inflated taxa in a hurdle step analysis
        a. Hurdle taxa = those taxa with 5% of more individuals with zero counts
        b. Zero samples set to “0”, non-zero samples set to “1” to for a test of presence vs absence
    - Lude Franke meta analysis pipeline.

The following Rscript is provided to assist in transformation of the taxa abundance and alpha-diversity data.
```
fgfp_gwas_cookbook_help.R 
```
It runs the following steps:
   - Identify the taxa with >=5% of samples with zero-values. 
     - While the above step would be the typical procedure we would like for you define the same taxa that we did as the hurdle taxa. These taxa are listed in the file “TaxaNamesHurdleOnly.txt”.
     - Define these as hurdle taxa.
     - create a binary variable for these taxa.
       - 0 = 0
       - !0 =1
     - create zero-truncated abundance data for these taxa.
       - all zero values are turned to NA and all non-zero values are treated as normal.
   - All remaining abundance data (both non-hurdle and hurdle-truncated abundance data) are rank normal transformed and log2 transformed
     - Rank normal transformation: completed with the rntransfrom() function from the GenABEL package
     - Log2 transformation: all 0 values are turned into NA, effectively performing a truncation of zero values just as done in the hurdle step.
   - Residualize RNT and LOG2 data
     - Fit a linear model in R, or your favorite programing language, setting the abundance values as the response and the following variables as explanatory variables. 
       - Top 10 PCs
       - Sex, age
       - Study specific batch variables, such as plate, processing date, etc..
     - Extract the residuals of the fitted model
     - Take care to return NAs for those samples that were already NAs. This is done to maintain order and structure in the data file for the association analysis.

#### 3 - Perform the association in SNPTEST.
  - Construct the .sample file for SNPTEST that includes all of the 
    - Sample ids
    - Covariates
    - LOG2 residual data
    - RNT residual data
    - Hurdle Binary data
  - Extract Variants to test
    - The ~25K variant positions to test are in the file “Sites_2_Test.txt”
    - This file contains HRCv1.1 SNPs and 1000G indels, as such you may not identify all of the variants in your data. This is OK! Proceed with all of those that you can identify.
      - We find 23,098 significant variants in the HRCv1.1 data.
    - Extract and prepare genotype data with qctools. Below is some code that could be used to do that:

```
for i in ~/YOURDIRECTORYPATH/*.bgen; do echo ${i}; qctool -g ${i} -incl-rsids SNPids2test.txt -ofiletype gen -og sigsites_chr${i# YOURDIRECTORYPATH /data_chr}.gen; done 

wc -l sigsites_chr*.gen  # 23,098 sites identified

## Concatenate all the sigsites together
for i in sigsites_chr*.gen; do echo ${i}; cat ${i} >> all_sigsites2text.gen; done
wc -l all_sigsites2text.gen # 23,098 sites 

## remove chromosome data
rm sigsites_chr*.gen
```
  - Perform association analysis for each taxa transformation

```
## RUN SNPTEST for Hurdle Binary Analysis
VAR1=C_Alphaproteobacteria_HB 
snptest_v2.5 -data all_sigsites2text.gen TransResCovar_MYSTUDY.sample -o $VAR1.out -log $VAR1.log -frequentist 1 -method score -cov_names PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 batchvar1 batchvar2 sex age -pheno $VAR1 -use_raw_phenotypes

##RUN SNPTEST for linear RNT and LOG analysis
VAR1=G_Blautia_RNT_Residuals
snptest_v2.5 -data all_sigsites2text.gen TransResCovar_MYSTUDY.sample -o $VAR1.out -log $VAR1.log -frequentist 1 -method score -pheno $VAR1 -use_raw_phenotypes
```

## (iv) Manhattan Plots for the association tests.
### Description
Running the following scripts will create a series of Manhattan plots, displaying the significance of each SNP associated with the microbial taxa tested. Plots for all the models run are produced and results can be merged by taxonomical group (eg. Proteobacteria) or by taxonomical level (eg. genera).

### Methods
#### 1 - Pull out the required information from tables with all the results.
For creating the Manhattan plots, only the SNP names ("SNP"), chromosomes ("CHR"), base pair coordinates ("BP") and p-values ("P") are needed (one p-value for each statistical model used to perform the MGWAS analysis).

```
# In a directory containing all the results for each taxa:
ls > all_taxa.txt
mkdir tables_mhplots && cd tables_mhplots/
bash produce_tables.sh
```

#### 2 - x
#### 3 - y
#### 4 - z

## (v) Phylogenetic trees in iTOL
### Description
Produce a phylogenetic tree of the taxa under study with overlaying piecharts containing information on loci counts (using different methods) and branches coloured according to heritability traits.

### Method:
#### 1 - Construct a phylogenetic tree using all the taxa in the study
Get all the taxa to include (genera, families, classes…) and load them in phyloT (http://phylot.biobyte.de/). This program generates phylogenetic trees based on the NCBI taxonomy. Choose options “expanded” and “polytomy yes”, otherwise random bifurcated structure is generated. 

**note**: Some names have to be corrected to NCBI taxonomy.

For example:
```
## From this list of taxa:
Bacteroidetes, Barnesiella, Bifidobacterium, Clostridia, Dorea, Prevotella

## phyloT generates this tree:
(((Prevotella,Barnesiella)Bacteroidales,(Dorea,Bifidobacterium)Terrabacteria_group)Bacteria);


                        /-Prevotella
           /Bacteroidales
          |             \-Barnesiella
-- /Bacteria
          |                   /-Dorea
           \Terrabacteria_group
                              \-Bifidobacterium


```

Some clades may present many polytomies - These can be resolved constructing trees based on the 16S rRNA gene sequences from RDP database.

Genus and next higher taxonomic levels are represented in dichotomies when they should be inclusive in a single branch.

#### 2 - Retrieve information on the loci counts
These can be obtained from the SigLocusCount table.
```
for tax in `cat all_taxa.txt`
do
   if grep -q $tax"_HurdleTruncRNT_Residuals" SNPbasedHeritability_SigLocusCount_v1_BOLTLMM.txt
      then grep $tax"_HurdleTruncRNT_Residuals" SNPbasedHeritability_SigLocusCount_v1_BOLTLMM.txt | cut -f 1,4
   else 
      echo -e $tax" \t0"
   fi
done > counts.txt
```
*to be updated with new tables*
#### 3 - Retrieve information on the heritability
*to be updated with new tables*
```
command
```
#### 4 - Upload data into iTOL and export tree
Run mgwas_itol.py script parsing the tree file and the tables for loci counts and heritability.

```
python3 mgwas_itol.py all_taxa.tree counts.txt heritability.txt
```
