#############################################
##  FGFP sample inference from amplicon data 
##  by Raul Tito & Rodrigo Bacigalupe
##  date: March 28th 2018
#############################################
### Pipeline created from https://benjjneb.github.io/dada2/tutorial.html
### More information available at https://benjjneb.github.io/dada2/index.html

### This pipeline requires demultiplexed files: Two fastq files per sample (*R1 and *R2), they can be *.gz

### Load required modules
# module load R/3.4.2   ### very important to use this version

### Load necessary libraries and set seed
library(dada2); packageVersion("dada2")
set.seed(12345)

### Every new MiSeq or Hiseq run is quality checked for number of reads and presence of unused barcodes to estimate 
### levels of carryover from prevous runs or contamination events. As part of this initial analysis demultipled files 
### per each sample are generated (Using LotuS).

### Set directory and files
path <- "~/../../data/DIndex_MiSeq_Runs/*year*/***RUN_NAME***/demultiplexed/" # CHANGE to the directory containing your demultiplexed R1 and R2 fastq files
fns <- list.files(path)
#fns

####################
### load your data
####################
fastqs <- fns[grepl(".fq.gz$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
### make sure that R1 is for forward read and R2 for reverse

fnFs <- fastqs[grepl(".1.fq.gz", fastqs)] ## Just the forward read files
fnRs <- fastqs[grepl(".2.fq.gz", fastqs)] ## Just the reverse read files
## Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, ".1.fq.gz"), `[`, 1) ## check if it is 1 or 2

## Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)


###########################################
## Examine quality profiles of F & R reads
###########################################
pdf("plotQualityProfile.pdf", onefile=T)
plotQualityProfile(fnFs[1:2]) ## remove 20 plus 10 (primers and first odd bases)
plotQualityProfile(fnRs[1:2]) ## remove 20 plus 10 (primers and first odd bases)
dev.off()

##################################
## Perform filtering and trimming
##################################
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

## Filter  the forward and reverse reads:
## Important to remove primers and low quality regions
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(130,200), ## Different settings tried,these are good for current primer constructs for MiSeq and HiSeq
                     trimLeft=c(30, 30),
                     maxN=0, maxEE=c(2,2), truncQ=11, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) #
head(out)

## Examine quality profiles of filtered reads
pdf("plotQualityProfile.filt.pdf", onefile=T)
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])
dev.off()


#########################
## Learn the Error Rates
#########################
## Learn forward error rates
errF <- learnErrors(filtFs, nread=1e6, multithread=TRUE) ## variable but this is the minimum number of reads
## Learn reverse error rates
errR <- learnErrors(filtRs, nread=1e6, multithread=TRUE) ## variable but this is the minimum number of reads
## Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

## Plot estimated error as sanity check 
pdf("plotErrors_F.pdf", onefile=T)
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf("plotErrors_R.pdf", onefile=T)
plotErrors(errR, nominalQ=TRUE)
dev.off()


#########################
## Perform dereplication
#########################
## Dereplicate the filtered fastq files
derepRs <- derepFastq(filtRs, verbose=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


####################
## Sample Inference
####################
## Apply the core sequence-variant inference algorithm to the dereplicated data
## Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

## Inspect the dada-class object returned by dada
dadaFs[[1]]

## Merge the denoised forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
## Inspect the merger data.frame from the first sample
head(mergers[[1]])


############################
## Construct sequence table 
############################
seqtab <- makeSequenceTable(mergers)
## Get dimensions
dim(seqtab)

## Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


###################
## Remove chimeras
###################
## Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


####################################
## Track reads through the pipeline
####################################
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)


###################
## Assign taxonomy
###################
taxHS <- assignTaxonomy(seqtab.nochim, "~/opt/dada2_db/rdp_train_set_16.fa.gz", multithread=TRUE) ## CHANGE to directory and pertinent database
colnames(taxHS) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
unname(head(taxHS))
unname(tail(taxHS))


#######################
## Write data to files
#######################
write.table(track, file = "track.tsv", quote=FALSE)
write.table(seqtab.nochim, file = "sequence_table_SV.tsv", quote=FALSE)
write.table(taxHS, file = "taxa_SV.tsv", quote=FALSE)
