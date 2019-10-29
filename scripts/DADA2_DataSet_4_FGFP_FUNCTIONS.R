###########################################
## DADA2 GWAS Data Preparation
## **** FUNCTIONS ****
##  by David Hughes
##  date: June 28th 2018
###########################################

### load necessary libraries
library(GenABEL)
library(vegan)
library(MASS)

###########################
## Microbiome SumStats
##  FUNCTION
##
############################
TaxaSumStatFun = function(taxa){
  d = as.numeric( na.omit( unlist(taxa) ) )
  ## proportion zero
  NZero = sum(d == 0)
  zeroProp = NZero/length(d)
  if(zeroProp > 0.95){
    out = c( zeroProp ,rep(NA,11) )
    names(out) = c("PropZero", "mean","sd","raw_Wstat",
                   "raw_Wpval", "Hurdle","trunc_mean","trunc_sd", 
                   "log2_Wstat", "log2_Wpval","rnt_Wstat","rnt_Wpval")
    } else{
  ##
      ## mean and sd
      msd = c(mean(d), sd(d))
      ## untransformed data normality test
      untrans_sh = shapiro.test(d)
      ## hurdle taxa ?
      hurdle = ifelse(zeroProp >= 0.05 & zeroProp <= 0.95, TRUE, FALSE)
      ##transformed normality test
      if(hurdle == FALSE){
        rnt_trans_sh = shapiro.test(rntransform(d))
        ## for the log transformation we must remove zeros
        f = d; f[f == 0] = NA;
        log_trans_sh = shapiro.test(log2(f))
        trunc_msd = c(NA, NA)
      } else {
        ## if the taxa is a hurdle taxa we will truncate the zero data 
        d[d == 0] = NA
        trunc_msd = c(mean(na.omit( d ) ), sd( na.omit(d) ))
        log_trans_sh = shapiro.test(log2(d))
        rnt_trans_sh = shapiro.test(rntransform(d))
      }
      ## DATA OUT
      out = c(zeroProp, msd, untrans_sh$statistic, untrans_sh$p.value, 
               hurdle, trunc_msd, log_trans_sh$statistic, 
              log_trans_sh$p.value, rnt_trans_sh$statistic, 
              rnt_trans_sh$p.value )
      names(out) = c("PropZero", "mean","sd","raw_Wstat",
                     "raw_Wpval","Hurdle","trunc_mean","trunc_sd", 
                     "log2_Wstat", "log2_Wpval","rnt_Wstat","rnt_Wpval")
      ###
    }
  return(unlist(out))
  }
  


###########################
## Load taxa to be included
##  in GWAS
##
############################
taxa2analyze = read.table("analyzable_taxa.txt", h = F, as.is = TRUE)[,1]
hurdletaxa = read.table("hurdletaxa.txt", h = F, as.is = TRUE)[,1]
## remove the alpha diversity
w = grep("Div_", taxa2analyze)
if(length(w)>0){
  taxa2analyze = taxa2analyze[-w]
  }
  


