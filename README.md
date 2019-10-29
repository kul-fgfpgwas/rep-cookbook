# rep-cookbook
Replication Cookbook for microbiome GWAS

### Microbiome processing: Sample inference from amplicon data
Perform sample inference from high-throughput amplicon data using [dada2](https://github.com/benjjneb/dada2).

The dada2 pipeline performs the following steps:
1. Quality control: Perform filtering and trimming of reads.
2. Perform dereplication: Group amplicon reads with same sequence into unique sequences.
3. Sample inference: Using the core sequence-variant inference algorithm to the dereplicated data.
4. Remove chimeras: Identify sequences that are exact bimeras.
5. Taxonomical classification.

Perform these analysis in R using the following script:
```
Seq_dada2.R
```
Transform the dada2 output into count tables and perform rarefaction using the following script:
```
dada2_to_taxtables.R
```

### Genotype QC and imputation

1. FGFP samples were processed on two different arrays, Human Core Exome v1.0 (576 samples) and Human Core Exome v1.1 (2112 samples), resulting in 545, 535 markers in merged raw data set. 
2. Allele calling was performed by GenomeStudio.
3. SNP QC: cross check with sex information, remove Unmapped (chr0) Variants, update to RS ids and remove duplicates; remove variants with >5% missing, remove monomorphic alleles, remove HWE p<1e-5 and remove AT/CG ambiguous sites. 
4. Sample QC: remove sample with >5% SNPs missing, remove samples with heterozygosity <3SD of mean value, remove samples with relatedness >0.025. 
5. Merge with 1KG samples and remove those with outstanding population structure. 
6. Imputation with HRCv1.1 server and EUR population as reference. 
7. Post imputation QC: remove SNPs that are monomorphic, INFO < 0.5 & MAF < 1% .

The scratch code in the following file will help to prepare the genotype data and prepare the files for the next step:
```
set_and_run_GWAS.txt
```

### Association tests
#### 1 - Select overlapping DADA2 between FGFP and Kiel cohorts.
There are a total of 92 taxa and 3 alpha diversity measures. DADA2 list can be found here:
```
analyzable_taxa.txt
```

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
       - Study specific batch variables (exttype, extractyear, aliquotyear, aliquotby, plate)
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

### Manhattan Plots for the association tests.

The following scripts will help to create the Manhattan plots displaying the significance of each SNP associated with the microbial taxa tested.

Pull out the required information from tables with all the results. For creating the Manhattan plots, only the SNP names ("SNP"), chromosomes ("CHR"), base pair coordinates ("BP") and p-values ("P") are needed (one p-value for each statistical model used to perform the MGWAS analysis).

```
# In a directory containing all the results for each taxa:
cat data_chr*out | grep -v "#" | grep -v "rs" | awk '{if ($19 >= 0.01 && $9 >= 0.3 && $22 >= 0.3) print
$2,$21}' | sed 's/:/\t/g' | sed 's/_.* /\t/g' > rnt.pvalues
```
Then use the following scripts to help you obtaining the manhattan plots:
```
mhplots.R
```

### Phylogenetic trees 
Produce a phylogenetic tree of the taxa under study with overlaying piecharts containing information on loci counts (using different methods) and branches coloured according to heritability traits.

#### Construct a phylogenetic tree using all the taxa in the study
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

Some clades may present many polytomies. These can be resolved constructing trees based on the 16S rRNA gene sequences from RDP database.

#### Retrieve information on the loci counts and heritability

These can be obtained from the SigLocusCount table.
```
for tax in `cat all_taxa.txt`
do
   if grep -q $tax"_HurdleTruncRNT_Residuals" SNPbasedHeritability_SigLocusCount_v1_BOLTLMM.txt
      then grep $tax"_HurdleTruncRNT_Residuals" SNPbasedHeritability_SigLocusCount_v1_BOLTLMM.txt | cut -f 1,4
   else 
      echo -e $tax" \t0"
   fi
done > loci_counts.txt
```

#### Upload data into iTOL and export tree
Load your .tree and .txt files with loci counts and heritability to https://itol.embl.de/ and format your tree [as desired](https://itol.embl.de/help.cgi)
