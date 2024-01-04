######################
######################
###### BAYESCAN ######
######################
######################

# script extracted from: 
https://rpubs.com/lbenestan/outlier


# 1- Downloads package

library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(qvalue)

# 2- prepare gest or bsc (bayescan file) to run into Bayescan program

vcf <- read.vcfR("RNAseq_mako_555loci_121samples.vcf")
pop_map <- read.table("pop_mako_RNAseq.txt", header=TRUE, stringsAsFactors = TRUE)
genind <- vcfR2genind(vcf)
genind@pop <- pop_map$STRATA
hierfstat <- genind2hierfstat(genind)
write.bayescan(hierfstat)


#3- run bayescan in a terminal

## Identifying outliers from bayescan analysis using R

bayescan=read.table("RNAseq_fst.txt") 
SNPb=read.table("id-555snp.txt",header=FALSE)
bayescan=cbind(SNPb, bayescan) 
colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST") 
write.table(bayescan, "mako_bayescan_results.txt", quote=FALSE, sep="\t", row.names=FALSE) 


attach(bayescan)
class(bayescan$Q_VALUE)
bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE) 
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 


bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4)) 
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4)) 
bayescan$ALPHA <- (round(bayescan$ALPHA, 4)) 
bayescan$FST <- (round(bayescan$FST, 6))

bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.05,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing")) 
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION)

positive <- bayescan[bayescan$SELECTION=="diversifying",] 
neutral <- bayescan[bayescan$SELECTION=="neutral",] 
balancing <- bayescan[bayescan$SELECTION=="balancing",]

# check the number of SNPS belonging to each category (neutral, balacing, selection)
xtabs(data=bayescan, ~SELECTION) 

#save the results of SNPS under selection and neutral in tables
write.table(neutral, "neutral.txt", row.names=F, quote=F)
write.table(balancing, "balancing.txt", row.names=F, quote=F) 
write.table(positive, "positive.txt", row.names=F, quote=F)
write.table(bayescan, "mako_loci_category_neutral_selection.txt", quote=FALSE, sep="\t", row.names=FALSE) 


# Transformation Log of the Q value in order to create the ggplot graph.
range(bayescan$Q_VALUE)
bayescan$LOG10_Q <- -log10(bayescan$Q_VALUE) 


# Make a nice graph w. ggplot2
x_title="Log(q-value)" 
y_title="Fst" 

tiff('mako_RNAseq_bayescan.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
ggplot(bayescan,aes(x=LOG10_Q,y=FST)) +
  geom_point(aes(fill=SELECTION), pch=21, size=2)+ 
  scale_fill_manual(name="Selection",values=c("white","red","orange"))+ 
  labs(x=x_title)+ 
  labs(y=y_title)+   
  theme_classic()
dev.off()


# the source file should be downloaded from Bayescan program
source("plot_R.r")
plot_bayescan("RNAseq_fst.txt")

