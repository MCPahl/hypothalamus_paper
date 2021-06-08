#Visualization marker genes
library(ggplot2)
dir="/mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/HypothalamusESC/DE"
setwd(dir)

###Graph a few important genes as line graphs, w/ standard deviation plotted as error bars

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

#Summarize data, calculate average and sd for tpm for each gene in each condition
tpm_avsd<- summarize_data(tpm,"tpm",c("class","gene_name"))

#Convert to factor so ggplot will put the samples in the right order
tpm_avsd$class <- factor(
	tpm_avsd$class,
	levels = c("HumanESC",
			   "HypothalamicProg",
			   "HypothalamicNeurons")
)

#Import list of genes
genelist <- read.delim("markerGenes.txt",header=FALSE,stringsAsFactors=FALSE)
markerGenes <- tpm_avsd[tpm_avsd$gene_name %in% genelist[,2],]
markerGenes$gene_name <- factor(markerGenes$gene_name, levels = genelist[,2])


#Generate line graphs for a few "marker genes"
col = "darkblue"
pdf("HypothalamusESC_rnaseq_HypothalamusMarkerGeneTPM.pdf")
ggplot(markerGenes, aes(x=class, y=tpm, group=gene_name, stat = "identity")) + 
    geom_errorbar(aes(ymin=tpm-sd, ymax=tpm+sd), width=.1, color=col) +
    geom_line(color=col) + geom_point(color=col)+
    theme_minimal() +
    theme(panel.grid.minor.x=element_blank(),
    	panel.grid.major.x=element_blank(),
    	panel.border = element_rect(colour = "black", fill=NA, size=1))+
    scale_x_discrete(labels=c("HumanESC" = "ESC", "HypothalamicProg" = "HP",
                              "HypothalamicNeurons" = "HN"))+
    labs(x = "", y = "Expression (TPM)" )+
    facet_wrap(gene_name ~ .,ncol = 3, scales = "free")
dev.off()
####