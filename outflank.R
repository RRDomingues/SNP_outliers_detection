#####################################
#####################################
############## OUTFLANK #############
#####################################
#####################################

# websites consult
# http://rstudio-pubs-static.s3.amazonaws.com/305384_9aee1c1046394fb9bd8e449453d72847.html
# https://tomjenkins.netlify.app/2020/09/21/r-popgen-getting-started/

# load packages
install.packages ("devtools")
library (devtools)
library(qvalue)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("qvalue")
library(vcfR)
library(ggplot2)
install.packages("patchwork")
library(patchwork)
library(adegenet)
library(hierfstat)
remotes::install_github("whitlock/OutFLANK")
library(OutFLANK)
install_github("green-striped-gecko/dartR")
library(dartR)
install.packages("VennDiagram")
library(VennDiagram)
library(viridis)

# read inputfile

vcf <- read.vcfR("RNAseq_mako_555loci_121samples.vcf")
pop_map <- read.table("pop_mako_RNAseq.txt", header=TRUE, stringsAsFactors = TRUE)
genind <- vcfR2genind(vcf)
genind@pop <- pop_map$STRATA
hierfstat <- genind2hierfstat(genind)
write.bayescan(hierfstat)


# Run OutFLANK using dartR wrapper script
outflank <- gl.outflank(genind,
                       withOutliers = TRUE,
                       NoCorr = TRUE,
                       Hmin = 0.05,
                       binwidth = 0.005,
                       Zoom = FALSE,
                       RightZoomFraction = 0.05,
                       titletext = NULL)


# Extract OutFLANK results
outflank.df = outflank$outflank$results


# Remove duplicated rows for each SNP locus
rowsToRemove = seq(1, nrow(outflank.df), by = 2)
outflank.df = outflank.df[-rowsToRemove, ]


# Print number of outliers (TRUE)
outflank.df$OutlierFlag %>% summary


# Extract outlier IDs
outlier_indexes = which(outflank.df$OutlierFlag == TRUE)
outlierID = locNames(genind)[outlier_indexes]
outlierID
write.table(outlierID, "outlier_RNAseq_outflank.txt")

#Convert Fsts <0 to zero
outflank.df$FST[outflank.df$FST < 0] = 0 


# Italic labels
fstlab = expression(italic("F")[ST])
hetlab = expression(italic("H")[e])


# Plot He versus Fst
tiff('mako_RNAseq_outflank.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
p.out.fk <- ggplot(data = outflank.df)+
  geom_point(aes(x = He, y = FST, colour = OutlierFlag))+
  scale_colour_manual(values = c("black","red"),
                      labels = c("Neutral SNP","Outlier SNP"))+
  xlab(hetlab)+
  ylab(fstlab)+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  size = 15,
                                  face = "bold"))
print(p.out.fk)
dev.off()

### Results
# outliers loci found: [1] "cp_152_TRINITY_DN115646_c10_g2_i2"  "cp_236_TRINITY_DN107499_c4_g3_i1"  
# [3] "cp_254_TRINITY_DN110114_c5_g3_i18"  "cp_472_TRINITY_DN112300_c2_g1_i2"  
# [5] "cp_504_TRINITY_DN119002_c5_g1_i1"   "cp_591_TRINITY_DN119887_c10_g1_i3" 
# [7] "cp_1060_TRINITY_DN119732_c2_g1_i1"  "cp_1192_TRINITY_DN119643_c1_g1_i8" 
# [9] "cp_1204_TRINITY_DN112142_c1_g4_i2"  "cp_1289_TRINITY_DN104894_c8_g1_i4" 
# [11] "cp_1308_TRINITY_DN119305_c6_g3_i1"  "cp_1755_TRINITY_DN117397_c5_g1_i6" 
# [13] "cp_2427_TRINITY_DN119509_c0_g2_i1"  "cp_2662_TRINITY_DN112988_c11_g1_i1"
# [15] "cp_2755_TRINITY_DN119835_c1_g1_i1"  "cp_3044_TRINITY_DN106563_c0_g1_i4" 
# [17] "cp_3398_TRINITY_DN118176_c5_g1_i4"  "cp_3574_TRINITY_DN110652_c3_g2_i3" 
# [19] "cp_5990_TRINITY_DN112009_c2_g2_i1" 
