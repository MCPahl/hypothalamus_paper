library("GenomicRanges")
library("tidyverse")

setwd("/mnt/isilon/sfgi/pahlm/analyses/grant/atacSeq/HypothalamusESC/annotation/promoteromeOverlap/")
load("promoteromeOCR_interaction.RData")

##Get a count of OCRs in promoters per celltype
promoterOCRs.split <- split(promoterOCRs.df, promoterOCRs.df$celltype)
promoter_ocr.count <- unlist(lapply(promoterOCRs.split, function(x){
	length(unique(x$OCR_consensusID))
	}))

##Get a count of promoter interacting OCRs per celltype
OCR2captureC.split <- split(OCR2captureC.df,OCR2captureC.df$OCR_celltype)
OCR2captureC.count <- unlist(lapply(OCR2captureC.split, function(x){
	length(unique(x$OCR_consensusID))
	}))

#Total nonpromoter OCRs, happens to be stored in this object, so decided to get the value from it
nearest.gene$celltype <- sub("_.*","", nearest.gene$OCR_name)
nearest.gene.split <- split(nearest.gene,nearest.gene$celltype)
npOCR_total.count <- unlist(lapply(nearest.gene.split, function(x){
	length(unique(x$OCR_consensusID))
	}))

nonPromoter_noCaptureC = npOCR_total.count - OCR2captureC.count

OCR_counts <- data.frame(OCR_promoter = promoter_ocr.count, 
			OCR_distal_promoterInteracting = OCR2captureC.count, 
			OCR_distal_nonInteracting = npOCR_total.count)
OCR_counts$celltype <- row.names(OCR_counts)

OCR_counts <- as_tibble(OCR_counts)
OCR_counts = OCR_counts %>% gather(key=OCR_class,value=count,OCR_promoter,OCR_distal_promoterInteracting,OCR_distal_nonInteracting)

OCR_counts$celltype = factor(OCR_counts$celltype, levels=c("HumanESCells","HypothalamicProg","HypothalamicNeurons"))
pdf("OCR_class_counts.pdf")
ggplot(OCR_counts, aes(x=celltype,y=count,group=OCR_class))+
theme_minimal() +
  geom_bar(stat="identity",position="dodge",width=0.75, fill="white", aes(col=celltype))+
  theme(panel.grid.minor.x=element_blank(),
            panel.grid.major.x=element_blank())
dev.off()
