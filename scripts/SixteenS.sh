###############################
#### Miseq processing pipeline
#### Jun Wang
#### Version 5th November 2014
###############################
#### This pipeline deals with V4 region from Miseq, 250*2 runs
#### You can run individual steps or the whole script as a .sh and submit to any core

#### Specifiy paths
MISEQ <- ~/data/Miseq_bin

#### Make demultiplexing easier after getting all the fastq (if ever zipped)
gunzip *gz


####################
#### Demultiplexing
####################
### MID.list looks like: (remove the ##, and for VDP: use _ or nothing between VDP and the number)
### SAMPLE	Index	GROUP
### SAMPLE1	Index1	TEST
### SAMPLE2	Index2	TEST

### Use the fastq-multx function from ea-ultilities 3 minutes
$MISEQ/ea-utils.1.1.2-686/fastq-multx -l MID.list Index.fastq Run1.fastq Run2.fastq -o n/a -o %.R1.fastq -o %.R2.fastq

### BETTER GET ALL RUNS BEFORE CREATING ANYTHING !!! ###
### remove useless fastq files
rm unmatched.*

### Make merge shell script:
### A shell script merging all reads per sample
ls *R1.fastq | sed 's/\./ /g' |uniq |awk '{print "$MISEQ/FLASH-1.2.10/flash -o "$1" "$1".R1.fastq "$1".R2.fastq -M 250 -m 230 -t 8"}' > flash.sh
### Speed up using parallel 1.5 min
cat flash.sh | $MISEQ/parallel-20140622/src/parallel -j 16
### Remove useless by-products from this script
rm *hist*
rm *notCombine*


### Quality filter shell
### change to fasta in one-go
### filtering threshold: 90% nucleotides must have quality score >25 
ls |grep extend |sed 's/\./ /g' |awk '{print "$MISEQ/fastx_toolkit-0.0.14/src/fastq_quality_filter/fastq_quality_filter -i "$1".extendedFrags.fastq -p 90 -q 25 -Q33 |$MISEQ/fastx_toolkit-0.0.14/src/fastq_to_fasta/fastq_to_fasta -o "$1".fasta -Q33"}' > quality_filter.sh
### speed up with parallel, 2.8 minutes
cat quality_filter.sh |$MISEQ/parallel-20140622/src/parallel -j 16


### Remove chimera from quality filtered reads 
### de-chime 
ls |grep fasta |sed 's/\./ /g' |awk '{print "~/software/usearch/6.0.307/usearch -uchime_ref "$1".fasta -db $MISEQ/gold.fa -strand plus -nonchimeras "$1".good.fasta"}' > dechime.sh
### parallel run
cat dechime.sh |$MISEQ/parallel-20140622/src/parallel -j 4


### Re-shuffle the sequence of each sample and take 10,000 reads out, give a new-id (10s) 
### Replace ID with artificial one, to mask the run information (important for submitting it to other database)
ls |grep good.fasta |sed 's/\./ /g'  |awk '{print "perl $MISEQ/subsample.pl "$1".good.fasta 10000 > "$1".sub10000.fasta"}' > subset.sh
### parallel run (10s) 
cat subset.sh | $MISEQ/parallel-20140622/src/parallel -j 16


### Classify 10,000 reads per sample (in 10 minutes)
ls |grep sub10000.fasta |sed 's/\./ /g'|awk '{print "java -jar $MISEQ/rdp_classifier_2.7/dist/classifier.jar -o "$1".tax -f fixrank "$1".sub10000.fasta"}' > classify.sh
### parallel run (in 10 mins) 
cat classify.sh | $MISEQ/parallel-20140622/src/parallel -j 16


### cat all tax file and produce taxonomical census
cat *tax > ALL.TAX

### make clean taxonomical table, 3 mins
perl $MISEQ/taxa_census_new.pl ALL.TAX

### Save space by zip fastq again
gzip *.fastq


### Make de-novo OTUs and create table
### usearch 6, small_mem, 20 min to 2 hours
### probably no upper limit for total sequences--haven't tested
### only to cross-check to genera results: so far genera is better in large population studies
### directly use perl to create an OTU table
cat *sub10000.fasta > ALL.sub10000.fasta
~/software/usearch/6.0.307/usearch -cluster_smallmem ALL.sub10000.fasta -centroids ALL.OTU.fasta -uc ALL.OTU.uc -id 0.97 -usersort
perl $MISEQ/uc_census.pl ALL.OTU.uc >OTU.tab.matrix

### Now can go on with stats with genus table and OTU table
