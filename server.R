library(shiny)
library(shinyjs)
library(DT)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(reshape2)
library(plotly)
library(Rsubread)

#setwd('/Users/valdezkm/Documents/Hackathon')
# load('./viSRA_practice_data.RData')

#### Make blast database ####
system(paste0('/home/ubuntu/ncbi-magicblast-1.3.0/bin/makeblastdb -in GCF_000001635.26_GRCm38.p6_genomic.fna -dbtype nucl -parse_seqids -out mouseGenome -title "mouseGenome "'))


shinyServer(function(input, output) {
  observeEvent(
    input$SRAbutton,
    isolate({ 
      write.csv(expression, file='normalized_data.csv')
      expression
    })
  )
  observeEvent(
    input$SRAbutton, {
    isolate({
      counts = reactive( 
        {
          SRA1 = input$SRRcode1
          SRA2 = input$SRRcode2
          
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
          GTF = '/home/ubuntu/Mus_musculus.GRCm38.91.gtf'
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
        
          dat <- counts()
          listGenes = read.delim(input$listOfGenes$datapath, sep = '\n', header = F)
          for (i in 1:length(listGenes[,1])) {
            geneName = listGenes[i,1]
            Genes = as.data.frame(dat[rownames(dat) %in% geneName,])
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
          v
        })
  
        output$sra=DT::renderDataTable(DT::datatable(
          {
            SRA1 = input$SRRcode1
            SRA2 = input$SRRcode2
            exp = counts()
            colnames(exp) = c(SRA1,SRA2)
            exp
          })
        )
        output$dotPlot <- renderPlotly({
          input$dotPlotButton
          isolate({
            dat <- as.data.frame(counts())
            d <- dat %>%
              mutate(geneID = rownames(dat)) %>%
              filter(geneID == toupper(input$geneName)) %>%
              melt() %>%
              ggplot(aes(x = variable, y = value)) +
              geom_point() +
              theme_bw() +
              xlab("SRA ID") +
              ylab("Expression (TPM)") +
              ggtitle(paste(input$geneName))
            d <- plotly_build(d)
            d$elementId <- NULL
            print(d)
          })
        })
        output$inc = renderUI(
          {
            tags$iframe(
              src='testing_fastqc.html',
              width='100%',
              height='800px')
          }
        )
    })
    })
})

