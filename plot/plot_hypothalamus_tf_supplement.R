library(tidyverse)

#load annotation
anno = read.delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/jasper2020/TF_motif.anno.txt",stringsAsFactors=FALSE)

#RNA-seq DE
DE.ESCvsHP = read.delim("/mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/HypothalamusESC/DE/tables/sigDE/ESCvsHP.sigDE_edgeR.txt",stringsAsFactors=F)
DE.HPvsHN = read.delim("/mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/HypothalamusESC/DE/tables/sigDE/HPvsHN.sigDE_edgeR.txt",stringsAsFactors=F)
load("/mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/HypothalamusESC/DE/tpm_data.Rdata")
tpm.m = tpm %>% group_by(gene_id, gene_name, class) %>% summarize(tpm = mean(log2tpm))
tpm.m$gene_id = gsub("\\..*","",tpm.m$gene_id)

#Transcription factor
tf = read.csv("bifet_enriched.table.csv", stringsAsFactors=F)
tf$FDR = NULL

tf.wide = tf %>% 
	as.tbl() %>%
	spread(key = celltype, value = p.bifet) %>% 
	mutate(ESC_rank = -log10(ESC)) %>%
	mutate(HP_rank = -log10(HP)) %>%
	mutate(HN_rank = -log10(HN)) %>%
	mutate(motif_id = gsub("^.*_","", TF))

x= left_join(tf.wide, anno, by="motif_id")
x.tpm = left_join(x, tpm.m, by= "gene_id")


pdf("ESC_HP_rank_scatterplot.pdf", useDingbats=FALSE)
ggplot(tf.wide, aes(x= ESC_rank, y= HP_rank))+
theme_minimal()+
geom_point()
dev.off()

pdf("HP_HN_rank_scatterplot.pdf", useDingbats=FALSE)
ggplot(tf.wide, aes(x=HN_rank , y= HP_rank))+
theme_minimal()+
geom_point()
dev.off()

 

 #gene sets 
#HP HN (shared top)
out = unique(x.tpm[x.tpm$HP_rank>10 & x.tpm$HN_rank>10,]$TF.y)
out = out[is.na(out)==F]
shared = tpm.m[tpm.m$gene_name %in% toupper(out),]
shared$class= factor(shared$class, levels= c("HumanESC", "HypothalamicProg", "HypothalamicNeurons"))
shared = shared %>% group_by(gene_name) %>% mutate(tmp.scaled = scale(tpm))

pdf("HP_HN_shared_tfs_heatmap.pdf", useDingbats=FALSE)
ggplot(shared, aes(x=class, y= gene_name, fill = tpm))+
	geom_tile(colour = "white") + 
	theme_minimal()+
	scale_fill_gradient(low = "white", high = "orange")
dev.off()


 #gene sets 
#HN more (shared top)
out = unique(x.tpm[x.tpm$HP_rank<10 & x.tpm$HN_rank>10,]$TF.y)
out = out[is.na(out)==F]
shared = tpm.m[tpm.m$gene_name %in% toupper(out),]
shared$class= factor(shared$class, levels= c("HumanESC", "HypothalamicProg", "HypothalamicNeurons"))
shared = shared %>% group_by(gene_name) %>% mutate(tmp.scaled = scale(tpm))

pdf("HN_enriched_notHP_tfs_heatmap.pdf", useDingbats=FALSE)
ggplot(shared, aes(x=class, y= gene_name, fill = tpm))+
	geom_tile(colour = "white") + 
	theme_minimal()+
	scale_fill_gradient(low = "white", high = "orange")
dev.off()

 #gene sets 
#HP more (shared top)
out = unique(x.tpm[x.tpm$HP_rank>10 & x.tpm$HN_rank<10,]$TF.y)
out = out[is.na(out)==F]
out[10:11] = c("ZIC1","ZIC2")
shared = tpm.m[tpm.m$gene_name %in% toupper(out),]
shared$class= factor(shared$class, levels= c("HumanESC", "HypothalamicProg", "HypothalamicNeurons"))
shared = shared %>% group_by(gene_name) %>% mutate(tmp.scaled = scale(tpm))

pdf("HP_enriched_notHP_tfs_heatmap.pdf", useDingbats=FALSE)
ggplot(shared, aes(x=class, y= gene_name, fill = tpm))+
	geom_tile(colour = "white") + 
	theme_minimal()+
	scale_fill_gradient(low = "white", high = "orange")
dev.off()


 #gene sets 
#HP more than ESC (shared top)
out = unique(x.tpm[x.tpm$HP_rank>2.5 & x.tpm$ESC_rank<10,]$TF.y)
out = out[is.na(out)==F]
shared = tpm.m[tpm.m$gene_name %in% toupper(out),]
shared$class= factor(shared$class, levels= c("HumanESC", "HypothalamicProg", "HypothalamicNeurons"))
shared = shared %>% group_by(gene_name) %>% mutate(tmp.scaled = scale(tpm))

pdf("HP_enriched_notESC_tfs_heatmap.pdf", useDingbats=FALSE)
ggplot(shared, aes(x=class, y= gene_name, fill = tpm))+
	geom_tile(colour = "white") + 
	theme_minimal()+
	scale_fill_gradient(low = "white", high = "orange")
dev.off()

