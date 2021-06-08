
###Analyze clusters for functional enrichment.
#library('org.Hs.eg.db')
load("/mnt/isilon/sfgi/pahlm/annotationFiles/msigdb_v7.0_GMTs/c5.bp.v7.0.symbols.Rdata")

HG.test = function(grp,i,anno){
  universe = unique(unlist(anno))
  gset = names(grp)[grp == i]
  gset = gset[gset %in% universe]
  out.hget = data.frame()
  
  tmp = lapply(anno,function(j){
      p = length(j [j %in% gset])
      m = length(j)
      n = length(universe) - m
      k = length(gset)
      hg.pval = phyper(p,m,n,k,lower.tail=FALSE)
    })
    out.hget = do.call(list,tmp)
}

go.enriched <-lapply(1:6,
  function(i){
  enrichment.hg<-HG.test(grp,i,c5.bp)
  enrichment.hg.adj<-p.adjust(enrichment.hg,method="BH")
  #enrichment.hg.adj[enrichment.hg.adj<0.05]
  enrichment.hg.adj.df<-as.data.frame(cbind(-log10(5e-67+enrichment.hg.adj[rank(enrichment.hg.adj)<11])))
  enrichment.hg.adj.df$GO.term = row.names(enrichment.hg.adj.df)
  row.names(enrichment.hg.adj.df) = NULL
  names(enrichment.hg.adj.df)[1] <- "Score"
  enrichment.hg.adj.df
})
names(go.enriched)<-paste("Group",1:6,sep="_")

#Dress the GO.term names up a little bit
rename<-lapply(go.enriched,function(group){
  z<-gsub("_"," ", sub("GO_","",group$GO.term))
  #nchar(z)
  z<- gsub("ENDOPLASMIC RETICULUM","ER",z)
  z<- gsub("NONSENSE MEDIATED DECAY","NMD",z)
  spaces.need = 71 - nchar(z)
  spaces <- lapply(spaces.need,function(x){
    paste(rep(" ", x),collapse="")
    })
  z<-paste(spaces,z)
  z
})

#Swap names
for(i in seq_along(rename)){
  go.enriched[[i]]$GO.term <- rename[[i]]
}

#Order by enrichment
go.enriched = lapply(go.enriched, function(go_group){
  go_group$GO.term = factor(go_group$GO.term, levels= go_group$GO.term[order(go_group$Score,decreasing = TRUE)])
  go_group
  })

#Build bargraphs for top 10 categories
go.plot.list<- lapply(go.enriched,function(group){
  ggplot(data=group, aes(x=GO.term, y=Score))+
    geom_bar(stat="identity",position="dodge",width=0.5,alpha=0.75)+
    scale_y_continuous('Length (mm)', limit=c(0,70))+
    theme_minimal() +
    theme(legend.position = "none")+
    theme(axis.text=element_text(size=6,
            family="sans"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    theme(panel.grid.minor.x=element_blank(),
            panel.grid.major.x=element_blank())+
    theme(plot.margin=unit(c(0,0,0,0),"cm"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    geom_abline(slope=0, intercept=1.30103,  col = "red",lty=2)
})
go.plot.grob<-lapply(go.plot.list,ggplotGrob)

pdf("GOterm_test_BH.pdf")
grid.arrange(grobs=go.plot.grob,
            ncol=6)
dev.off()
