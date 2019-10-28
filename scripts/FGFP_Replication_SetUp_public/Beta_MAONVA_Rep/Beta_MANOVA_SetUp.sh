#################################
## BETA-DIVERSITY MANOVA
##  HELPER SCRIPTS
##
## by: David Hughes
## date: July 17th 2018
#################################


#######################
## Step 1:
##  Make Necessary Dosage
##  Files
#######################
## chromosome(s) to screen (chr -> 2,3,4,6,9,16)
array=(02 03 04 05 09 16)
for i in "${array[@]}"
	do 
	CHR=data_chr${i}.vcf
	echo ${CHR}
	##
	qctool -g ${CHR} -filetype vcf -vcf-genotype-field GP  -incl-rsids beta_sites_2_test.txt -og ${CHR}.dosage -ofiletype dosage
done

#######################
## Step 2:
##  cat dosage files
##  together
#######################
touch FOCUS_betasites2test.dosage
for i in *.dosage; do cat ${i} >> FOCUS_betasites2test.dosage; done



##################################
##  RUN Beta GWAS !
##  Date: July 2nd, By David Hughes
##################################

## Load "phenodata" Phenotype/Covariate Data Matrix
	## mine was already in a Rdata file
	## as the element "phenodata"
	## where everything was str() checked to be
	## factors or numeric
load("MYSTUDY_covariate_matrix.RData")
	
	## ALTERNATIVE: read a flat text file of "phenodata" if need be
	# phenodata = read.table("", header = TRUE, sep = "\t")

## LOAD DATA
	#args = commandArgs(trailingOnly=TRUE)
## file that holds the genotype 
	#genotypefile = args[1]
genotypefile = "FOCUS_betasites2test.dosage"

## Load Genetic Dta
gendata = read.table(file = genotypefile, header = T, sep = " ")
## ****GENOTYPE START AT COLUMN 7 (BE SURE TO CHECK THIS IS TRUE !!!)*****
gencolumns = 7:ncol(gendata)
sampleids = colnames(gendata)[gencolumns]

#################################################
### FUNCTION TO PERFORM MANOVA
#################################################
MANOVAfun = function(genotype, covardata){
  genotype = unlist(genotype)
  ### MANOVA 1
  ### ****** MAKE THE MODEL YOUR OWN !!!! ******
  Man0 <- manova( betadiv ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                    PC7 + PC8 + PC9 + PC10 + snpsex + age + 
                    genotype , data = covardata)
  #Wstat = apply(Man0$residuals, 2, function(x){shapiro.test(x)$statistic})
  ss =  summary(Man0)[[4]]
  Man0_out = ss[nrow(ss)-1, c(2,6)]
  ## data out
  out = c(Man0_out)
  names( out ) = c( "fit0_Pillai", "fit0_pvalue" )
  return(out)
  }

#################################################
### END OF FUNCTION TO PERFORM MANOVA
#################################################

######################
### EXECUTE ANALYSIS
######################
manovaout = t( apply(gendata[, gencolumns], 1, function(x){
  out = MANOVAfun(genotype = x, covardata = phenodata )
  return(out)
  }) )
## add snpids as rownames
rownames(manovaout) = gendata[,3]

######################
### Write data 2 file
######################
n = paste0( genotypefile, "_betadiv_manova.out")
write.table(manovaout, file = n, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)




