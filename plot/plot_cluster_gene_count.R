############Gene types ###########################
load("/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/source/gene_tss.RData")
#Loads gene_len, gene_type, gencodeV19
grp.table <- as.tibble(data.frame(group=grp,gene=names(grp)))
gene_type %>% select(gene,type_detail) -> gene_type.g
gene_type.g <- unique(gene_type.g)

grp.gene_types <- left_join(grp.table,gene_type.g)
grp.gene_types %>% group_by(group) %>% summarise(group_count=n()) -> grp.group.summary
grp.gene_types %>% group_by(group,type_detail) %>% summarise(count=n()) -> grp.gene_types.summary
grp.gene_types.summary <- left_join(grp.gene_types.summary, grp.group.summary)
grp.gene_types.summary <- grp.gene_types.summary[
  grp.gene_types.summary$type_detail %in% c("lincRNA","protein_coding"),
]

#Get differences
grp.gene_types.summary %>% 
  group_by(group) %>% 
  summarise(group_count=mean(group_count),diff_count=sum(count)) %>% 
  mutate(count=group_count-diff_count) %>% 
  select(group,count,group_count) %>%
  add_column(type_detail="Other") -> grp.gene_types.others

grp.gene_types.summary<- bind_rows(grp.gene_types.summary, grp.gene_types.others)


grp.gene_types.summary <- grp.gene_types.summary %>% mutate(ratio = count/group_count*100) 


grp.gene_types.summary$group <- factor(
  grp.gene_types.summary$group,
  levels = c(2,6,4,3,5,1)
)

grp.gene_types.summary$type_detail <- factor(
  grp.gene_types.summary$type_detail,
  levels = c("protein_coding", "lincRNA", "Other")
)


grp.table %>% group_by(group) %>% summarise(count=n()) -> grp.counts

grp.counts$group <- factor(
  grp.counts$group,
  levels = c(2,6,4,3,5,1)
)

pdf("Cluster_Counts.pdf")
ggplot(grp.counts, aes(x=group,y=count,group=group))+
theme_minimal() +
  geom_bar(stat="identity",position="dodge",width=0.75, aes(fill=group))+
  theme(panel.grid.minor.x=element_blank(),
            panel.grid.major.x=element_blank())+
  scale_fill_manual("legend",values=c("1" = "#999999", "2" = "#56B4E9", "3" = "#009E73", "4" = "#CC79A7", "5" = "#F0E442", "6" = "#E69F00" ))
dev.off()

### Plot 
pdf("Cluster_GeneType.pdf")
ggplot(grp.gene_types.summary, aes(x=type_detail,y=ratio,group=group))+
theme_minimal() +
  geom_bar(stat="identity",position="dodge",width=0.75, aes(fill=group))+
  theme(panel.grid.minor.x=element_blank(),
            panel.grid.major.x=element_blank())+
  scale_fill_manual("legend",values=c("1" = "#999999", "2" = "#56B4E9", "3" = "#009E73", "4" = "#CC79A7", "5" = "#F0E442", "6" = "#E69F00" ))
dev.off()
