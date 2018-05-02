# rep-cookbook
Replication Cookbook for FGFP GWAS

## (i) Microbiome processing: Sample inference from amplicon data
### Description
Perform sample inference from high-throughput amplicon data using the SixteenS pipeline and [dada2](https://github.com/benjjneb/dada2).

**SixteenS pipeline performs the following tasks:**
1. Merging reads: Forward and reverse reads are merged with the FLASH software (v1.2).
2. Reads QC: quality filtering was subsequently performed with the fastx toolkit (http://hannonlab.cshl.edu/fastx_toolkit/), excluding those sequences with >5% nucleotides of quality score <30.
3. Chimera removal: chimeras in sequences were removed using UCHIME (v6.0) with gold.fa as reference. 
4. Subsampling: randomly selecting 10,000 reads for each sample.
5. Taxonomical classification and compositional matrices for each taxonomical level were carried out using RDP classifier with RDP14, where classifications with low confidence at the genus level (<0.8) were organized in an arbitrary taxon of "unclassified_family".

**dada2 pipeline performs the following steps:**
1. Merging reads: Forward and reverse reads are merged with the FLASH software (v1.2).

### Methods:
Run analysis using the SixteenS pipeline (in bash):
```
bash SixteenS.sh
```

Run analysis using dada2 (in R):
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

## (iii) Sample inference from amplicon data
