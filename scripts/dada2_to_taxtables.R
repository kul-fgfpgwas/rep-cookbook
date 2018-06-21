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


###################
# Import datasets #
###################

#Import 16S pipeline dataset
#Import genus CMM data used for GWAS analysis
taxa_16s_CMM <- read.table("../GENERA_CMM_renamed.tsv", header=TRUE, row.names = 1)

#Remove 0 after "VDP." in samplenames
rownames(taxa_16s_CMM) <- gsub("VDP.0", "VDP.", rownames(taxa_16s_CMM))

#Import mapping table for VDP IDs
labID <- read.table("../LabID_VDP.txt", header=TRUE)

#Change names
ind <- match(row.names(taxa_16s_CMM), labID$VDP_.ID)
row.names(taxa_16s_CMM) <- labID$LabID[ind]


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
psdada2 <- dada2_to_phyloseq("../FlinHS010203_dada2.otu_table.txt", "../FlinHS010203_dada2.tax_table.txt",sep="\t")

#Remove samples with low read counts (below 10000)
num_reads <- sample_sums(psdada2)[sample_sums(psdada2) < 10000]
Sample_names <- names(num_reads)
psdada2.above10000 = prune_samples(names(which(sample_sums(psdada2) > 10000)), psdada2)

#Get names of control samples (samples with Bl, EXBL, NC, ... in names)
out.samples<-names(sample_sums(psdada2.above10000))
out.samples<- c(
  out.samples[ grep("BL",out.samples) ],
  out.samples[ grep("Bl",out.samples) ],
  out.samples[ grep("EXBL",out.samples) ],
  out.samples[ grep("NC",out.samples) ],
  out.samples[ grep("NP",out.samples) ],
  out.samples[ grep("RS",out.samples) ],
  out.samples[ grep("Zymmo",out.samples)],
  out.samples[ grep("_pl_",out.samples)],
  out.samples[ grep("PC",out.samples) ] 
)

#Replace manually some sample names
colnames(otu_table(psdada2.above10000)) <- sub("VDP_138_Rep.r01","VDP_138",colnames(otu_table(psdada2.above10000)))
colnames(otu_table(psdada2.above10000)) <- sub("VDP_47_Rep.r01","VDP_47",colnames(otu_table(psdada2.above10000)))
colnames(otu_table(psdada2.above10000)) <- sub("VDP_700_Rep_HiSeq01","VDP_700",colnames(otu_table(psdada2.above10000)))
colnames(otu_table(psdada2.above10000)) <- sub("VDP_2895_Rep","VDP_2895",colnames(otu_table(psdada2.above10000)))

#Remove control samples
in.samples <- names(table(c(out.samples, colnames(otu_table(psdada2.above10000))))[table(c(out.samples, colnames(otu_table(psdada2.above10000)))) == 1])
psdada2.above10000.noctrl = prune_samples(in.samples, psdada2.above10000)

#Remove ".r01" extension in dada2 sample names
colnames(otu_table(psdada2.above10000.noctrl)) <- sub(".r0.*","",colnames(otu_table(psdada2.above10000.noctrl)))

#Keep only samples that are also in 16S pipeline CMM for comparison
J16s.samples <- row.names(taxa_16s_CMM)
psdada2.above10000.noctrl.filt = prune_samples(J16s.samples, psdada2.above10000.noctrl)

#Summarize to the genus taxa level using tax_glom
psdada2.gen.glom <- tax_glom(psdada2.above10000.noctrl.filt, taxrank = "Genus") #This step may take a while

#Perform rarefaction at 10000
psdada2.gen.glom.rar = rarefy_even_depth(psdada2.gen.glom, replace=FALSE, sample.size = 10000)

#Write genus rarefied table
tax.table.gen.rar <- otu_table(psdada2.gen.glom.rar)
rownames(tax.table.gen.rar) <- tax_table(psdada2.gen.glom.rar)[rownames(tax.table.gen.rar),6]
write.table(tax.table.gen.rar, "Genus_dada2_Rar.tsv" ,col.names=T,row.names = T,quote=FALSE,sep = "\t") 

#Write family rarefied table
tax.table.fam.rar <- otu_table(psdada2.gen.glom.rar)
rownames(tax.table.fam.rar) <- tax_table(psdada2.gen.glom.rar)[rownames(tax.table.fam.rar),5]
#Sum up all the rows with same taxa names
tax.table.fam.rar <- aggregate(as.data.frame(tax.table.fam.rar), list(row.names(tax.table.fam.rar)), sum) 
write.table(tax.table.fam.rar, "Family_dada2_Rar.tsv" ,col.names=T,row.names = ,quote=FALSE,sep = "\t")

#Write order rarefied table
tax.table.ord.rar <- otu_table(psdada2.gen.glom.rar)
rownames(tax.table.ord.rar) <- tax_table(psdada2.gen.glom.rar)[rownames(tax.table.ord.rar),4]
#Sum up all the rows with same taxa names
tax.table.ord.rar <- aggregate(as.data.frame(tax.table.ord.rar), list(row.names(tax.table.ord.rar)), sum) 
write.table(tax.table.ord.rar, "Order_dada2_Rar.tsv" ,col.names=T,row.names = T,quote=FALSE,sep = "\t")

#Write class rarefied table
tax.table.cla.rar <- otu_table(psdada2.gen.glom.rar)
rownames(tax.table.cla.rar) <- tax_table(psdada2.gen.glom.rar)[rownames(tax.table.cla.rar),3]
#Sum up all the rows with same taxa names
tax.table.cla.rar <- aggregate(as.data.frame(tax.table.cla.rar), list(row.names(tax.table.cla.rar)), sum) 
write.table(tax.table.cla.rar, "Class_dada2_Rar.tsv" ,col.names=T,row.names = T,quote=FALSE,sep = "\t")

#Write phylum rarefied table
tax.table.phy.rar <- otu_table(psdada2.gen.glom.rar)
rownames(tax.table.phy.rar) <- tax_table(psdada2.gen.glom.rar)[rownames(tax.table.phy.rar),2]
#Sum up all the rows with same taxa names
tax.table.phy.rar <- aggregate(as.data.frame(tax.table.phy.rar), list(row.names(tax.table.phy.rar)), sum) 
write.table(tax.table.phy.rar, "Phylum_dada2_Rar.tsv" ,col.names=T,row.names = T,quote=FALSE,sep = "\t")
