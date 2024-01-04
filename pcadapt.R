######################
######################
###### PCADAT ########
######################
######################

#https://bcm-uga.github.io/pcadapt/articles/pcadapt.html

#install packages
library(pcadapt)
library(vcfR)
library(qvalue)

# transform VCF to BED format using PLINK program
# plink --vcf RNAseq_mako_555loci_121samples.vcf --make-bed --allow-extra-chr 0 --out mako_555loci_121samples

# read BED input file
data <- read.pcadapt("mako_555loci_121samples.bed", type = "bed")

#Choosing the number k of principal components
x <- pcadapt(input = data, K = 20)

#Scree plot
plot(x, option = "screeplot")

plot(x, option = "screeplot", K = 10)


#Score plot

# With integers
poplist.int <- c(rep(1, 4), rep(2, 17), rep(3, 16), rep(4, 4), rep(5, 31), rep (6, 15), rep (7, 34))

poplist.names <- c(rep("PTG", 4), rep("EAO", 17), rep("CV", 16), rep("BRZ", 4), rep("SAO", 31), rep ("SAF", 15), rep ("IO", 34))

print(poplist.int)

print(poplist.names)


plot(x, option = "scores", pop = poplist.int)

plot(x, option = "scores", pop = poplist.names)

plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)


## Compute the statistic test based on PCA

x <- pcadapt(data, K = 2)
summary(x)

#Graphic tools

#Manhattan plot
plot(x , option = "manhattan")

# Q-Q plot
plot(x, option = "qqplot", threshold = 0.1)

# Histograms of the test statistic and of the p-values
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

plot(x, option = "stat.distribution")

# Choosing a cutoff for outlier detection

#q-values
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)
length(outliers)

#Benjamini-Hochberg Procedure
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)

# Bonferroni correction
padj1 <- p.adjust(x1$pvalues,method="bonferroni")
alpha <- 0.1
outliers1 <- which(padj < alpha)
length(outliers)
outliers 



