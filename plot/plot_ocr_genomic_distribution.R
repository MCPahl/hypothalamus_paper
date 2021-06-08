##########
#ALL OCRS#
##########

library(GenomicRanges)
library(ggplot2)
library(tidyverse)
#Read my datafiles 
myATACdir = c("/mnt/isilon/sfgi/pahlm/analyses/grant/atacSeq/HypothalamusESC/comparison/quantitative/")
files =  list.files(myATACdir)
myBedFiles = files[grep("mergedOCRs_filtered_FPKMabove1.bed",files)]
myBedPaths = paste(myATACdir, myBedFiles, sep="/")
myOCRs = lapply(myBedPaths, function(myBedPath){
	read.delim(myBedPath, stringsAsFactors=FALSE, header=FALSE,sep=" ")
})

#Sample label
labels = c("ESC","HN","HP")
for(i in seq_along(labels)){
	myOCRs[[i]]$label = c(labels[i])
}

#Convert my OCRs to single dataframe, then convert to GenomicRanges object
myOCRs = do.call("rbind", myOCRs)
myOCRs.gr = GRanges(
	seqnames= myOCRs$V1,
	ranges = IRanges(myOCRs$V2,myOCRs$V3),
	OCR_ID = myOCRs$V4,
	strand = c("*"),
	label = myOCRs$label
)

#Read annotation (from Chun)
dir="/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/"
files= paste(dir,c("gencode.v19.annotation.transcript_only.TSS.bed",
	"genecode_v19_3UTR_comp.clean.bed",
	"genecode_v19_5UTR_comp.clean.bed",
	"genecode_v19_CDS_comp.clean.bed",
	"genecode_v19_intron_comp.clean.bed"),sep="")
beds <- lapply(files, function(file){
	read.delim(file, stringsAsFactors=FALSE,header=FALSE)
	})

#Give each feature name and combine to one dataframe
names(beds) = c("TSS", "gene_3UTR", "gene_5UTR", "gene_CDS","gene_intron")
beds = do.call("rbind",beds)
beds$annotation_type = gsub("\\..*","",row.names(beds))
row.names(beds)=NULL

#Change TSS locations to our lab's definition of promoter (2kb downstream of TSS)
beds$V2 = ifelse(beds$annotation_type=="TSS", ifelse(beds$V6=="+", beds$V2-1500, beds$V2), beds$V2)
beds$V3 = ifelse(beds$annotation_type=="TSS", ifelse(beds$V6=="+", beds$V3+500, beds$V3), beds$V3)

beds$V2 = ifelse(beds$annotation_type=="TSS", ifelse(beds$V6=="-", beds$V2-500, beds$V2), beds$V2)
beds$V3 = ifelse(beds$annotation_type=="TSS", ifelse(beds$V6=="-", beds$V3+1500, beds$V3), beds$V3)
beds$annotation_type[beds$annotation_type=="TSS"] = "Promoter"

#Convert annotation to GenomicRange object
beds.gr = GRanges(
	seqnames=beds$V1,
	ranges = IRanges(beds$V2,beds$V3),
	names = beds$V4,
	geneID = gsub("_.*","", gsub("\\+.*","",beds$V4)),
	value = beds$V5,
	strand = beds$V6,
	annotation_type=beds$annotation_type
)

beds.gr[beds.gr$annotation_type=="gene_intron"& grepl("_intron_1_",beds.gr$names),]$annotation_type = "gene_first_intron"

#Find overlap with features
overlap_index = as.data.frame(findOverlaps(myOCRs.gr, beds.gr))

OCR.feature_mapped = myOCRs.gr[overlap_index[,1]]
feature.OCR_mapped = beds.gr[overlap_index[,2]]
OCR2feature = data.frame(OCR.feature_mapped, feature.OCR_mapped)
feature_hier= c("Promoter", "gene_5UTR" ,"gene_3UTR", "gene_CDS", "gene_first_intron", "gene_intron")


OCR2feature.l = split(OCR2feature,OCR2feature$OCR_ID)
OCR_topfeature = lapply(OCR2feature.l, function(OCR2feature){
	ifelse(feature_hier[1] %in% unique(OCR2feature$annotation_type), feature_hier[1],
		ifelse(feature_hier[2] %in% unique(OCR2feature$annotation_type), feature_hier[2],
			ifelse(feature_hier[3] %in% unique(OCR2feature$annotation_type), feature_hier[3],
				ifelse(feature_hier[4] %in% unique(OCR2feature$annotation_type), feature_hier[4],
					ifelse(feature_hier[5] %in% unique(OCR2feature$annotation_type), feature_hier[5],feature_hier[6])
					)
				)
			)
		)
})

gene_topFeature = unlist(OCR_topfeature)
feature_anno = data.frame(feature= gene_topFeature, OCR=names(gene_topFeature))

intergenic = unique(myOCRs.gr$OCR_ID)[unique(myOCRs.gr$OCR_ID) %in% feature_anno$OCR==FALSE]
intergenic = data.frame(feature="intergenic" , OCR=intergenic)
feature_anno = as_tibble(rbind(feature_anno, intergenic))

feature_counts = feature_anno  %>% group_by(feature) %>% tally()
feature_counts$feature= factor(feature_counts$feature, levels=c("Promoter", "gene_5UTR", "gene_CDS","gene_first_intron", "gene_intron", "gene_3UTR", "intergenic"))
feature_counts$percentage = 100*feature_counts$n/sum(feature_counts$n)

pdf("OCR_distribution.pdf")
ggplot(feature_counts, aes(x="", y=percentage,fill=feature))+
	geom_bar(width=1, stat="identity")+
	coord_polar("y", start=0)+
	theme_minimal()
dev.off()

##########
#PIR-OCRs#
##########


library(GenomicRanges)
library(ggplot2)
library(tidyverse)
#Read my datafiles 
myATACdir = c("/mnt/isilon/sfgi/pahlm/analyses/grant/atacSeq/HypothalamusESC/annotation/promoteromeOverlap/promoter_interacting_region_annotations/")
files =  list.files(myATACdir)
myBedFiles = files[grep("_mergedOCRsANDPromoters_filtered_FPKMabove1_noExtend.bed",files)]
myBedPaths = paste(myATACdir, myBedFiles, sep="/")
myOCRs = lapply(myBedPaths, function(myBedPath){
	read.delim(myBedPath, stringsAsFactors=FALSE, header=FALSE,sep="\t")
})

myOCRs <- lapply(myOCRs, function(myOCR){
myOCR[grep("PIR",myOCR$V6),]
})

#Sample label
labels = c("ESC","HN","HP")
for(i in seq_along(labels)){
	myOCRs[[i]]$label = c(labels[i])
}

#Convert my OCRs to single dataframe, then convert to GenomicRanges object
myOCRs = do.call("rbind", myOCRs)
myOCRs.gr = GRanges(
	seqnames= myOCRs$V1,
	ranges = IRanges(myOCRs$V2,myOCRs$V3),
	OCR_ID = myOCRs$V4,
	strand = c("*"),
	label = myOCRs$label
)

#Read annotation (from Chun)
dir="/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/"
files= paste(dir,c("gencode.v19.annotation.transcript_only.TSS.bed",
	"genecode_v19_3UTR_comp.clean.bed",
	"genecode_v19_5UTR_comp.clean.bed",
	"genecode_v19_CDS_comp.clean.bed",
	"genecode_v19_intron_comp.clean.bed"),sep="")
beds <- lapply(files, function(file){
	read.delim(file, stringsAsFactors=FALSE,header=FALSE)
	})

#Give each feature name and combine to one dataframe
names(beds) = c("TSS", "gene_3UTR", "gene_5UTR", "gene_CDS","gene_intron")
beds = do.call("rbind",beds)
beds$annotation_type = gsub("\\..*","",row.names(beds))
row.names(beds)=NULL

#Change TSS locations to our lab's definition of promoter (2kb downstream of TSS)
beds$V2 = ifelse(beds$annotation_type=="TSS", ifelse(beds$V6=="+", beds$V2-1500, beds$V2), beds$V2)
beds$V3 = ifelse(beds$annotation_type=="TSS", ifelse(beds$V6=="+", beds$V3+500, beds$V3), beds$V3)

beds$V2 = ifelse(beds$annotation_type=="TSS", ifelse(beds$V6=="-", beds$V2-500, beds$V2), beds$V2)
beds$V3 = ifelse(beds$annotation_type=="TSS", ifelse(beds$V6=="-", beds$V3+1500, beds$V3), beds$V3)
beds$annotation_type[beds$annotation_type=="TSS"] = "Promoter"

#Convert annotation to GenomicRange object
beds.gr = GRanges(
	seqnames=beds$V1,
	ranges = IRanges(beds$V2,beds$V3),
	names = beds$V4,
	geneID = gsub("_.*","", gsub("\\+.*","",beds$V4)),
	value = beds$V5,
	strand = beds$V6,
	annotation_type=beds$annotation_type
)

beds.gr[beds.gr$annotation_type=="gene_intron"& grepl("_intron_1_",beds.gr$names),]$annotation_type = "gene_first_intron"

#Find overlap with features
overlap_index = as.data.frame(findOverlaps(myOCRs.gr, beds.gr))

OCR.feature_mapped = myOCRs.gr[overlap_index[,1]]
feature.OCR_mapped = beds.gr[overlap_index[,2]]
OCR2feature = data.frame(OCR.feature_mapped, feature.OCR_mapped)
feature_hier= c("Promoter", "gene_5UTR" ,"gene_3UTR", "gene_CDS", "gene_first_intron", "gene_intron")


OCR2feature.l = split(OCR2feature,OCR2feature$OCR_ID)
OCR_topfeature = lapply(OCR2feature.l, function(OCR2feature){
	ifelse(feature_hier[1] %in% unique(OCR2feature$annotation_type), feature_hier[1],
		ifelse(feature_hier[2] %in% unique(OCR2feature$annotation_type), feature_hier[2],
			ifelse(feature_hier[3] %in% unique(OCR2feature$annotation_type), feature_hier[3],
				ifelse(feature_hier[4] %in% unique(OCR2feature$annotation_type), feature_hier[4],
					ifelse(feature_hier[5] %in% unique(OCR2feature$annotation_type), feature_hier[5],feature_hier[6])
					)
				)
			)
		)
})

gene_topFeature = unlist(OCR_topfeature)
feature_anno = data.frame(feature= gene_topFeature, OCR=names(gene_topFeature))

intergenic = unique(myOCRs.gr$OCR_ID)[unique(myOCRs.gr$OCR_ID) %in% feature_anno$OCR==FALSE]
intergenic = data.frame(feature="intergenic" , OCR=intergenic)
feature_anno = as_tibble(rbind(feature_anno, intergenic))

feature_counts = feature_anno  %>% group_by(feature) %>% tally()
feature_counts$feature= factor(feature_counts$feature, levels=c("Promoter", "gene_5UTR", "gene_CDS","gene_first_intron", "gene_intron", "gene_3UTR", "intergenic"))
feature_counts$percentage = 100*feature_counts$n/sum(feature_counts$n)

pdf("PIR_distribution.pdf")
ggplot(feature_counts, aes(x="", y=percentage,fill=feature))+
	geom_bar(width=1, stat="identity")+
	coord_polar("y", start=0)+
	theme_minimal()
dev.off()

