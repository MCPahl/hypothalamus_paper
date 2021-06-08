library(tidyverse)
library(GenomicRanges)

load("/mnt/isilon/sfgi/pahlm/analyses/grant/atacSeq/HypothalamusESC/annotation/promoteromeOverlap/promoteromeOCR_interaction.RData")
lept_pos = read.delim("GSE112125_lepr_pos_hg19_liftover.bed", stringsAsFactors=F, header=F)
lept_pos = lept_pos[,c(1:3)]
names(lept_pos) = c("chr", "start", "end")

HN = OCR2captureC.df[OCR2captureC.df$OCR_celltype=="HypothalamicNeurons",]
names(HN)[c(1,2,3)]=c("chr", "start","end")


index = as.data.frame(findOverlaps(GRanges(data.frame(HN)), GRanges(lept_pos)))

HN.lept.pos = data.frame(HN[index[,1],])

lept.pos.overlap = unique(data.frame(HN[index[,2],]))


HN_lept.pos_nb2b = HN.lept.pos[HN.lept.pos$prey_feature==".",]
HN_lept.pos_nb2b = HN_lept.pos_nb2b %>% mutate(dist = (bait_start+bait_end)/2 - (prey_start+prey_end)/2)

lepr_pos_hist_overlap = unique(gsub("N.*\\+","",unlist(strsplit(HN_lept.pos_nb2b$bait_feature,split="\\|"))))


load("/mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/HypothalamusESC/DE/tpm_data.Rdata")
tpm = tpm[tpm$class=="HypothalamicNeurons",]
tpm = tpm %>% group_by(gene_id, gene_name, class) %>% summarise(tpm = mean(tpm))

tpm = tpm[tpm$gene_name %in% lepr_pos_hist_overlap,]
#tpm_expressed = tpm[tpm$tpm > 1,]
tpm_expressed = tpm
HN_lep.pos.tpm = tpm_expressed[tpm_expressed$gene_name %in% lepr_pos_hist_overlap,]


HN_lep.pos.tpm$condition = "HN_lepR.pos"

pdf("lept_pos_H3k27ac_connected_gene.pdf", width =3)
ggplot(HN_lep.pos.tpm, aes(x= class, y= log2(tpm+1)))+
 geom_violin(fill="orange")+
 theme_minimal()+
 geom_hline(yintercept=log2(2))+
 theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()


nrow(unique(lept_pos))
nrow(lept.pos.overlap)

percentage = data.frame(class = c("overlap_ocr", "non_overlap"), percentage = c(28.7, 71.3))

pdf("H3k27ac_inersect_connected_gene.pdf", width =3)
ggplot(percentage, aes(fill=class, y=percentage, x=1)) + 
    geom_bar(position="stack", stat="identity")+
     theme_minimal()+
 theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()



