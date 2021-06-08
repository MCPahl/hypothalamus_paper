library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(ggplot2)
library(gridExtra)
library("Biostrings")
library("Rsamtools")
library("GenomicRanges")
library("tidyverse")

###################################
#Load clusters and RNA-seq tpm Data

#Load promoterOCRs (promoterOCRs.df), promoter interacting OCRs (OCR2captureC.df), and nearest gene mappings for each OCR (nearest.gene)
load("/mnt/isilon/sfgi/pahlm/analyses/grant/atacSeq/HypothalamusESC/annotation/promoteromeOverlap/promoteromeOCR_interaction.RData")

#Load ATAC-seq FPKM data
fpkm <- as_tibble(read.delim("/mnt/isilon/sfgi/pahlm/analyses/grant/atacSeq/HypothalamusESC/comparison/quantitative/fpkm_unfiltered.txt",
					sep=" ",
					header=TRUE,
					stringsAsFactors=FALSE))
#fix up fpkm
fpkm = fpkm[,-c(1)]
fpkm = fpkm %>% dplyr::rename("OCR_consensusID"=id) 
###################################
#Combine promoterOCRs with 
promoterOCRs.df
promoterOCRs_fpkm = left_join(promoterOCRs.df,fpkm)

#function to summarize data
summarize_data <- function(data, varname, groupnames){
  require(plyr)
  summarize_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum <- ddply(data, groupnames, .fun=summarize_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

fpkm = fpkm %>% gather(key=sample,value=fpkm, 
	"HumanESCells_rep1","HumanESCells_rep2", "HumanESCells_rep3",
	"HumanESCells_rep4","HypothalamicProg_rep1", "HypothalamicProg_rep2",  
 	"HypothalamicProg_rep3",    "HypothalamicProg_rep4","HypothalamicNeurons_rep1",
 	"HypothalamicNeurons_rep2","HypothalamicNeurons_rep3", "HypothalamicNeurons_rep4",
	"HypothalamicNeurons_rep5", "HypothalamicNeurons_rep6")

fpkm$class = sub("_rep.*","",fpkm$sample)

fpkm_avsd<- summarize_data(fpkm,"fpkm",c("class","OCR_consensusID"))

#Convert to factor so ggplot will put the samples in the right order
fpkm_avsd$class <- factor(
	fpkm_avsd$class,
	levels = c("HumanESCells",
			   "HypothalamicProg",
			   "HypothalamicNeurons")
)

OCR_list <- "OCR_107216"

NKX2.1_promoter <- fpkm_avsd[fpkm_avsd$OCR_consensusID  %in% OCR_list,]

#Generate line graphs for a few "marker genes"
col = "darkblue"
pdf("HypothalamusESC_ATAC-seq_NKX2-1_PromoterAccessibility.pdf")
ggplot(NKX2.1_promoter, aes(x=class, y=fpkm, group=OCR_consensusID, stat = "identity")) + 
    geom_errorbar(aes(ymin=fpkm-sd, ymax=fpkm+sd), width=.1, color=col) +
    geom_line(color=col) + geom_point(color=col)+
    theme_minimal() +
    theme(panel.grid.minor.x=element_blank(),
    	panel.grid.major.x=element_blank(),
    	panel.border = element_rect(colour = "black", fill=NA, size=1))+
    scale_x_discrete(labels=c("HumanESC" = "ESC", "HypothalamicProg" = "HP",
                              "HypothalamicNeurons" = "HN"))+
    labs(x = "", y = "Accessibility (FPKM)" )+
    facet_wrap(OCR_consensusID ~ .,ncol = 3, scales = "free")
dev.off()
