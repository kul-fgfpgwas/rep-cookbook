# Scripts for converting dada2 output into count tables for all taxa levels

#############################
# Install and load packages #
#############################

source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
biocLite('microbiome')
library(phyloseq)
library(dplyr)
library(microbiome)

#Set seed (for reproducibility)
set.seed(1234)

####################
# Define functions #
####################

#Function to add unclassified categories at genus level (got it from Jorge)
dada2_to_phyloseq = function(abundance.file,tax.file,sep=" "){
  library("phyloseq"); packageVersion("phyloseq")
  library("data.table")
  #system.time(data<-read.table(abundance.file, header=T, row.names=1, dec=".", sep=sep))
  system.time(data <- fread(abundance.file, header=T, fill = TRUE, dec=".", sep=sep) )
  data.colnames<-colnames(data)[1:length(colnames(data))-1]
  data <-data.frame(data, row.names=1)
  colnames(data) = data.colnames	
  
  data<-t(data)
  tax<-read.table(tax.file, header=T, row.names=1, dec=".", sep=sep)
  table(rownames(data) == rownames(tax))
  tax<-as.matrix(tax)
  
  for(i in rownames(tax)){
    temp.line<-tax[i,]
    temp.line[1]<-paste("K",temp.line[1],sep="_")
    temp.line[2]<-paste("P",temp.line[2],sep="_")
    temp.line[3]<-paste("C",temp.line[3],sep="_")
    temp.line[4]<-paste("O",temp.line[4],sep="_")
    temp.line[5]<-paste("F",temp.line[5],sep="_")
    temp.line[6]<-paste("G",temp.line[6],sep="_")
    #if( length(grep("T", is.na(temp.line))) > 0 ){ #### Looks for the NA (missing taxonomic labels)
    if(length(grep("_NA",  temp.line )) > 0 ){ #### Looks for the NA (missing taxonomic labels)
      na.ind<-grep("_NA",  temp.line )  ### Locates in the row where are such values
      na.lev<-min( na.ind )                   ### Locates the first NA taxonomic label
      last.tax<-temp.line[na.lev-1]           ### Takes the last-known taxonomic category
      ###last.tax<-as.matrix(last.lev)[1]         
      temp.line[ na.ind ] <- rep( paste( "unclassified",last.tax, sep="_" ),length(temp.line[ na.ind ]) ) ### Add an uc_LAST_KNOWN_TAX_LEVEL to all the unkwon taxonomic levels, ie: 
      tax[i,]<-temp.line     ### Change the line for the one whith the new taxonomic info
    } else{
      tax[i,]<-temp.line     ### Change the line for the one whith the new taxonomic info
    }
    
  }
  
  OTU = otu_table(data, taxa_are_rows = TRUE)
  TAX = tax_table(tax[rownames(OTU),])
  physeq = phyloseq(OTU, TAX)
  return(physeq)
  
}
#End of function

####################
# Perform analysis #
####################

#Load dada2 data using the dada2_to_phyloseq function (files produced by Raul, @ /raeslab/group/data/HiSeq3000plus_dada2/)
psdada2 <- dada2_to_phyloseq("../otu_table.txt", "../tax_table.txt",sep="\t")

#Remove samples with low read counts (below 10000)
num_reads <- sample_sums(psdada2)[sample_sums(psdada2) < 10000]
Sample_names <- names(num_reads)
psdada2.above10000 = prune_samples(names(which(sample_sums(psdada2) > 10000)), psdada2)

#Summarize to the genus taxa level using tax_glom
psdada2.gen.glom <- tax_glom(psdada2.above10000, taxrank = "Genus") #This step may take a while

#Perform rarefaction at 10000
psdada2.gen.glom.rar = rarefy_even_depth(psdada2.gen.glom, replace=FALSE, sample.size = 10000)

#Write taxa rarefied table (here genus)
tax.table.gen.rar <- otu_table(psdada2.gen.glom.rar)
rownames(tax.table.gen.rar) <- tax_table(psdada2.gen.glom.rar)[rownames(tax.table.gen.rar),6] #change this value for fam, order, class...
write.table(tax.table.gen.rar, "Genus_dada2_Rar.tsv" ,col.names=T,row.names = T,quote=FALSE,sep = "\t") 
