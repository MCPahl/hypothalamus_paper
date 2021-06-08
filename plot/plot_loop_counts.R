
library(ggplot2)
library(tidyverse)

files = list.files()
files = files[grepl("ibed", files)]

counts = lapply(files, function(file){
	frag_int = read.delim(file,stringsAsFactors=FALSE)
	count = nrow(unique(frag_int))
	data.frame(cell =  gsub("_chicagoloops_merged.ibed", "", file), count = count)
})
counts = do.call("rbind", counts)
counts$cell = factor(counts$cell, levels= c("ESC", "HP", "HN"))

write.table(counts, file="CaptureC_loop.counts.txt", quote=F, row.name=F, sep="\t")


pdf(file="chromatin_loops_captureC_count.pdf", useDingbats=FALSE)
ggplot(counts, aes(x = cell, y = count))+
 geom_point(color="darkblue",fill="darkblue", size=3)+
 theme_minimal()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top")+
    scale_y_continuous(breaks = seq(0, 550000, by= 100000))+
    expand_limits(y=0)+
    scale_x_discrete(limits = rev(levels(counts$cell)))+
    coord_flip()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()