listGenes = read.delim('ListGenes.txt', sep='\n', header = F)
for (i in 1:length(listGenes[,1])) {
  geneName = listGenes[i,1]
  Genes = v[rownames(v) %in% geneName,]
  Genes = data.frame(Genes)
  Genes$samples = rownames(Genes)
  jpeg(file = paste0(geneName,"_DotPlot.jpeg"))
  print(ggplot(Genes, aes(x=samples, y=Genes)) + 
          geom_dotplot(binwidth=.25, binaxis = 'y', stackdir = 'center', dotsize = 1) + 
          theme(panel.background = element_blank(), axis.text=element_text(size=12),
                axis.title=element_text(size=14), plot.title=element_text(size=20,hjust=0.5)) +
          labs(x='', y=paste0(geneName,' Expression'), size=5))
  dev.off()
}

ggplot(Genes, aes(x=samples, y=Genes)) + geom_violin() + geom_jitter(shape=16, position = position_jitter(0.2), size=1) +
  geom_boxplot(outlier.sie = -1, width = 0.25, fill = 'antiquewhite', color = 'antiquewhite', alpha = 0.5) + labs(title="Violin Plots", x = " ", y = "log2expression") + 
  theme(legend.position = 'none') + coord_flip()

heatmap.2(Genes, scale = 'row')