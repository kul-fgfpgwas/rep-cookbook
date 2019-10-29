##########################
# MGWAS - Manhattan Plots
##########################

# Create MHplots for specific analysis (RNT, LOG, HB)

# Import required libraries
require(data.table) # For efficient data handling
library(qqman) # Faor Manhattan plots
library(ggbio)
require(biovizBase)
require(GenomicRanges)
require(dplyr)
library(gdata)

#Import data and create tables
rnt <- read.table("RNT.pvalues", header=FALSE)
colnames(rnt) <- c("CHR", "BP", "P", "TAXA")
rnt$SNP <- paste(rnt$CHR, rnt$BP, sep=":")
rnt$source <- c("rnt")

log <- read.table("LOG.pvalues", header=FALSE)
colnames(log) <- c("CHR", "BP", "P", "TAXA")
log$SNP <- paste(log$CHR, log$BP, sep=":")

hb <- read.table("HB.pvalues", header=FALSE)
colnames(hb) <- c("CHR", "BP", "P", "TAXA")
hb$SNP <- paste(hb$CHR, hb$BP, sep=":")
hb$source <- c("hb")

#Merge the three datasets for estimating x-axis sizes
all <- rbind(rnt, hb)

#Create object to plot
gwas <- all %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(all, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) 
  
axisdf = gwas %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#Plot hits
ggplot(data=gwas) +
  
  # Show all points
  geom_point(data=subset(gwas, gwas$P < 0.000000000159), aes(x=BPcum, y=-log10(P), color=as.factor(TAXA)), alpha=1, size=1.8) +
  scale_color_manual(values=col) +
  geom_point(data=subset(gwas, gwas$P > 0.000000000159), aes(x=BPcum, y=-log10(P)), color="grey85", alpha=1, size=1.8) +
  geom_point(data=subset(gwas, gwas$P > 0.000000000159 & gwas$CHR %in% seq(1,22,2)), aes(x=BPcum, y=-log10(P)), color="grey75", alpha=1, size=1.8) +
   
  # custom X axis:
  scale_x_continuous(expand = c(0, 0), label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits=c(5,10), breaks=c(5,6,7,8,9,10)) +     # remove space between plot area and x axis
  
  geom_hline(yintercept=9.798603, linetype="dashed", size=1, color = "red") +
  geom_hline(yintercept=7.60206, linetype="dashed", size=1, color = "black") +
 #Axis legends
  xlab("Chromosome") + ylab("-log10(P value)") +

  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    axis.title.x = element_text(size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14)
  )

dev.off()

#meta-Analysis

metaanalysis <- read.table("snps_metaanalysis.txt", header=TRUE)
metaanalysis$SNP <- paste(metaanalysis$chr, metaanalysis$pos, sep=':')
metaanalysis$SNPTAX <- paste(metaanalysis$SNP, metaanalysis$tax, sep=':')
colnames(metaanalysis) <- c("TAXA", "METHOD", "CHR", "BP", "P", "SNP", "SNPTAX")

metaanalysis$BPcum <- metaanalysis$SNPTAX
metaanalysisMP <- merge(metaanalysis, ALL, by="SNPTAX")

svg("METAANALYSIS.svg", width = 14, height = 4)

ggplot(data=metaanalysis) +
  # Show all points
  geom_point(aes(x=BPcum, y=-log10(em_P_value), color=COLOR), alpha=0.8, size=2) +

  # custom X axis:
  scale_x_continuous(expand = c(0, 0), label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,14), breaks=c(2,4,6,8,10,12,14)) +     # remove space between plot area and x axis
  
  geom_hline(yintercept=9.798603, linetype="dashed", size=1, color = "black") +
  geom_hline(yintercept=7.60206, linetype="dashed", size=1, color = "gray") +
  geom_hline(yintercept=5, linetype="dashed", size=1, color = "red") +
  
  #Axis legends
  xlab("Chromosome") + ylab("-log10(P)") +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
    axis.title.x = element_text(size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14)
  ) 

dev.off()
