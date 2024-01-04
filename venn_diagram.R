#Website to edit ggvenndiagram
#https://www.datanovia.com/en/blog/beautiful-ggplot-venn-diagram-with-r/

#utilise the website above to know shared terms. 

# https://bioinfogp.cnb.csic.es/tools/venny/index.html

# Libraries
library(tidyverse)
library(VennDiagram)
library (ggVennDiagram)

# Load dataset from github
data <- read.table("venn_diagram_RNAseq_genome_scan.txt", header=TRUE) 

#just in case need to remove some variable
#to_remove <- c("_|[0-9]|\\.|function|^id|script|var|div|null|typeof|opts|if|^r$|undefined|false|loaded|true|settimeout|eval|else|artist")
#data <- data %>% filter(!grepl(to_remove, word)) %>% filter(!word %in% stopwords('fr')) %>% filter(!word %in% proust_stopwords()$word)

#select variables (in this case, program e loci)
x1 = list(
  data %>% filter(program=="Bayescan") %>% select(loci) %>% unlist(), 
  data %>% filter(program=="OutFLANK") %>% select(loci) %>% unlist(), 
  data %>% filter(program=="pcadapt") %>% select(loci) %>% unlist())
  
#plot veen diagram
ggVennDiagram (x1)

#put the names of variables
ggVennDiagram(x1[1:3], label_alpha = 0, category.names = c("Bayescan", "OutFLANK", "pcadapt"))

# remove to background and colored edge

tiff('venn_diagram_RNAseq_genome_scan.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
ggVennDiagram(
  x1, label_alpha = 0,
  category.names = c("Bayescan" , "OutFLANK" , "pcadapt"),
  edge_size = 0, label_color = "white",
) 
dev.off ()

# graph without percentage
tiff('venn_diagram_withoutercentage_RNAseq_genome_scan.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
ggVennDiagram(
  x1, label_alpha = 0,
  category.names = c("Bayescan" , "OutFLANK" , "pcadapt"),
  edge_size = 0, label_color = "white", label = "count",
) 
dev.off ()



