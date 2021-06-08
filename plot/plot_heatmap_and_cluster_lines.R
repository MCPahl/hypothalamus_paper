
##Heatmap of differentially expressed genes
library(ComplexHeatmap)
library(tidyverse)
#library(dendextend)
library()
library(ggplot2)
library(gridExtra)

#Get differentially expressed genes base
delist <- unique(sigDE_df$gene_name)
DE.tpmlist <- tpm[tpm$gene_name %in% delist,]
DE.tpmlist <- DE.tpmlist[,-c(5,6,7)]
DE.tpmlist %>% spread(sample,tpm) -> DE.tpmlist
DE.tpmlist <- as.data.frame(DE.tpmlist)
DE.tpmlist <- DE.tpmlist[duplicated(DE.tpmlist$gene_name)==FALSE,]
row.names(DE.tpmlist) <- DE.tpmlist[,2]
DE.tpmlist <- DE.tpmlist[,-c(1,2)]
DE.tpmlist <- as.matrix(DE.tpmlist)
DE.tpmlist <- DE.tpmlist[apply(DE.tpmlist,1,sum)>0,]

t(apply(DE.tpmlist,1,scale))->scaled.DE
colnames(DE.tpmlist)->colnames(scaled.DE)


#Get differentially expressed genes log2
DE.tpmlist.log2 <- tpm[tpm$gene_name %in% delist,]
DE.tpmlist.log2 <- DE.tpmlist.log2[,-c(4,6,7)]
DE.tpmlist.log2 %>% spread(sample,log2tpm) -> DE.tpmlist.log2
DE.tpmlist.log2 <- as.data.frame(DE.tpmlist.log2)
DE.tpmlist.log2 <- DE.tpmlist.log2[duplicated(DE.tpmlist.log2$gene_name)==FALSE,]
row.names(DE.tpmlist.log2) <- DE.tpmlist.log2[,2]
DE.tpmlist.log2 <- DE.tpmlist.log2[,-c(1,2)]
DE.tpmlist.log2 <- as.matrix(DE.tpmlist.log2)

#Maximum norminalization
#DE.tpmlist.max <- apply(DE.tpmlist,1,max)
#DE.tpmlist<- DE.tpmlist[!DE.tpmlist.max==0,]
#DE.tpmlist.max <- DE.tpmlist.max[!DE.tpmlist.max==0]
#DE.tpmlist.norm <- DE.tpmlist/DE.tpmlist.max

#Get the max for each column
DE.tpmlist.log.max <- apply(DE.tpmlist.log2,1,max)
DE.tpmlist.log.max <- as.matrix(DE.tpmlist.log.max)
DE.tpmlist.log.max <- DE.tpmlist.log.max[row.names(DE.tpmlist.log.max) %in% row.names(scaled.DE),]
DE.tpmlist.log.max <- as.matrix(DE.tpmlist.log.max)

#Clustering

#Hierarchal clustering, use pearson correlation and "agnes" clustering algorithm (two closest correlations grouped, iteratively)
clust<-hclust(as.dist(1-cor(t(scaled.DE), method="pearson")), method="ward.D2")

grp <- cutree(clust,k=6)


col_fun1 = colorRamp2(c(-2,-1,0,1,2),c("blue","skyblue","white","lightcoral","red"),space="RGB")
col_fun2 = colorRamp2(c(quantile(DE.tpmlist.log.max ,0), quantile(DE.tpmlist.log.max ,.9)),c("white","darkorange"),space="RGB")

#Prepare the first heatmap, scores normalized by maximum value
h_list = Heatmap(scaled.DE,
  show_row_names=FALSE,
  col=col_fun1,
  cluster_rows=clust,
  cluster_columns = FALSE,
  column_order=c(1,2,3,7,8,9,4,5,6),
  heatmap_legend_param = list(title = "Scaled Expression", col_fun=col_fun1,color_bar="continuous"),
  split = 6)+
Heatmap(grp, name = "clusters", show_row_names = FALSE, width = unit (5, "mm"),
    col = structure(names = c("1", "2", "3", "4","5","6"), 
      c("#999999","#56B4E9","#009E73" ,"#CC79A7","#F0E442" ,"#E69F00" )))+
Heatmap(DE.tpmlist.log.max,
  name = "max log2(TPM+1)",
  show_row_names=FALSE,
  col=col_fun2,
  heatmap_legend_param = list(title = "max log2(TPM+1)", col_fun=col_fun2,color_bar="continuous"))

pdf("heatmap.pdf")
draw(h_list)
dev.off()


###Plot the pattern of each cluster
DE.tpm.df <- as.data.frame(DE.tpmlist)
DE.tpm.df$gene_name <- row.names(DE.tpm.df)
row.names(DE.tpm.df) <- NULL
DE.tpm.df <- as.tibble(DE.tpm.df)

grp.df <- as.data.frame(grp)
grp.df$gene_name <- row.names(grp.df)
row.names(grp.df) <- NULL
grp.df <- as.tibble(grp.df)
row.names(DE.tpm.df) <- NULL

DE.tpm.df = left_join(DE.tpm.df,grp.df)

DE.tpm.df = DE.tpm.df %>% gather(key=sample,value=tpm, "HumanESC_N_rep1","HumanESC_N_rep2",           
 "HumanESC_N_rep3","HypothalamicNeurons_N_rep1",
 "HypothalamicNeurons_N_rep2","HypothalamicNeurons_N_rep3",
 "HypothalamicProg_N_rep1","HypothalamicProg_N_rep2",
 "HypothalamicProg_N_rep3")

grp.summary$class = sub("_N_.*","",DE.tpm.df$sample)

grp.summary <- DE.tpm.df %>% 
  group_by(grp,class) %>%
  summarise(
            quantile05 = log2(as.numeric(quantile(tpm,.05))),
            quantile25 = log2(as.numeric(quantile(tpm,.25))),
            median = log2(median(tpm)),
            quantile75 = log2(as.numeric(quantile(tpm,.75))),
            quantile95 = log2(as.numeric(quantile(tpm,.95)))
  )
  
  grp.summary<- transform(grp.summary,color.use = 
                                    ifelse(grp==1, "#999999",
                                      ifelse(grp==2,"#56B4E9",
                                        ifelse(grp==3, "#009E73",
                                          ifelse(grp==4, "#CC79A7",
                                            ifelse(grp==5, "#F0E442",
                                              "#E69F00")))))
  )

grp.summary$class <- factor(
  grp.summary$class,
  levels = c("HumanESC",
         "HypothalamicProg",
         "HypothalamicNeurons")
)
                                      
pdf("group_linegraph.pdf")
ggplot(grp.summary, aes(x=class, y=median, group=grp, col=color.use, stat = "identity")) + 
    geom_ribbon(aes(ymin=quantile25, ymax=quantile75,fill=color.use,linetype=NA,alpha=0.1))+
    geom_line(size=1)+
    geom_boxplot(width=0.25,
      mapping=aes(group = class,
      lower = quantile25,
      upper = quantile75,
      middle = median,
      ymin = quantile05,
      ymax = quantile95),
      stat="identity")+
    theme_minimal() +
    theme(panel.grid.minor.x=element_blank(),
      panel.grid.major.x=element_blank())+
    scale_x_discrete(labels=c("HumanESC" = "ESC", "HypothalamicProg" = "HP",
                              "HypothalamicNeurons" = "HN"))+
    labs(x = "", y = "Expression (TPM)" )+
    scale_color_manual(values = unique(grp.summary$color.use))+
    scale_fill_manual(values = unique(grp.summary$color.use))+
    facet_wrap(grp ~ .,ncol = 2)
dev.off()