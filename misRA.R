library(shiny)
library(shinyjs)
library(DT)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(reshape2)
library(plotly)
library(Rsubread)
library(limma)
library(edgeR)

#First argument should be the first SRR
#Second argument should be the second SRR
#Third argument should be text file of genes

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

# setwd('/Users/valdezkm/Documents/Hackathon')
#load('./viSRA_practice_data.RData')

#### Make blast database ####
system(paste0('/home/ubuntu/ncbi-magicblast-1.3.0/bin/makeblastdb -in GCF_000001635.26_GRCm38.p6_genomic.fna -dbtype nucl -parse_seqids -out mouseGenome -title "mouseGenome "'))

SRA1 = args[1]
SRA2 = args[2]
listGenes = read.delim(args[3])

#### Blast SRR numbers against reference genome ####
#SRR531311
#SRR531315
system(paste0('/home/ubuntu/ncbi-magicblast-1.3.0/bin/magicblast -sra ',SRA1,' -db /home/ubuntu/mouseGenome -no_unaligned -num_threads 8 -out sam_mouse1.sam'))
system(paste0('/home/ubuntu/ncbi-magicblast-1.3.0/bin/magicblast -sra ',SRA2,' -db /home/ubuntu/mouseGenome -no_unaligned -num_threads 8 -out sam_mouse2.sam'))
#system(paste0(magicDir,'/magicblast -sra ',SRA1,SRA2,' -db ',magicDir,'/GRCh38 -num_threads 8'))

#### SAM to BAM ####
system(paste0('samtools view -S -b sam_mouse1.sam > sam_mouse1.bam'))
system(paste0('samtools view -S -b sam_mouse2.sam > sam_mouse2.bam'))

#### Count reads with Rsubread ####
GTF = '/home/ubuntu/GCF_000001635.26_GRCm38.p6_genomic.gff'
system(paste0('/home/ubuntu/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene -a ',GTF,' -o counts1.txt sam_mouse1.bam'))
system(paste0('/home/ubuntu/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene -a ',GTF,' -o counts2.txt sam_mouse2.bam'))

#### Read in counts ####
system(paste0("awk '{if (NR!=1) {print}}' counts_gene_both.txt > counts_gene_both1.txt"))
counts_all = read.delim('counts_gene_both1.txt', sep="\t")
rownames(counts_all) = counts_all$Geneid
counts_all = subset(counts_all, select = c(sam_mouse1.bam,sam_mouse2.bam))

#### Filter Counts ####
keep <- rowSums(counts_all) >= 5
counts_all <- counts_all[keep,]

#### Normalize with TMM ####
dge = DGEList(counts_all)
dge = calcNormFactors(dge, method='TMM')
v <- voom(dge)

#### Dot Plots ####
#listGenes = read.delim('ListGenes.txt', sep='\n', header = F)
for (i in 1:length(listGenes[,1])) {
  geneName = listGenes[i,1]
  Genes = as.data.frame(v$E[rownames(v$E) %in% geneName,])
  colnames(Genes) = geneName
  jpeg(file = paste0(geneName,"_DotPlot.jpeg"))
  if (nrow(Genes)!=0) {
    Genes$Samples = c(SRA1, SRA2)
    print(ggplot(Genes, aes(x=Genes[,2], y=Genes[,1])) + 
            geom_dotplot(binwidth=.25, binaxis = 'y', stackdir = 'center', dotsize = 0.5) + 
            theme(panel.background = element_blank(), axis.text=element_text(size=12),
                  axis.title=element_text(size=14), plot.title=element_text(size=20,hjust=0.5)) +
            labs(x='', y=paste0(geneName,' Expression'), size=5))
  }
  dev.off()
}

#### Write out normalized counts ####
colnames(v$E) = c(SRA1,SRA2)
write.csv(v$E, file = 'normalized_counts.csv')

  

