##Plot venn diagram for groups

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(pathview))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(clusterProfiler))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="orthogroups.tsv' [default %default]",
              dest="orthogroups"),
  make_option(c("-o","--output"), type="character", default="all",
              help="output filename prefix [default %default]",
              dest="output_filename"),
  make_option(c("-a","--eggnog"), type="character", default=NULL,
              help="folder with 1v1 functional annotations [default %default]",
              dest="eggnog"),
  make_option(c("-b","--ancestral"), type="character", default=NULL,
              help="folder with ancestral state output [default %default]",
              dest="ancestral"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)  

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

##Manually add files
{
  setwd("/home/yichun/huilab/butterfly/r8s_lambda3/")
  opt$orthogroups="Orthogroups.tsv"
  opt$anno="/home/yichun/3enrichment"
}

orthogroups<-read.delim(opt$orthogroups,
                        sep = "\t", header = TRUE)

pvalue=0.05
p.table<-read.delim("Gamma_family_results.txt", header = TRUE)
p.table<-p.table[p.table$pvalue<pvalue,]
base.change.tab<-read.delim("Gamma_change.tab", header = FALSE)
names(base.change.tab)<-base.change.tab[1,]
base.change.tab<-base.change.tab[base.change.tab$FamilyID %in% p.table$X.FamilyID, ]
count.gl<-as.data.frame(names(base.change.tab[2:ncol(base.change.tab)]))
names(count.gl)[1]="cafe.label"
species.gain<-base.change.tab[as.numeric(base.change.tab[,2])>2,1]

ortho.butterfly<-orthogroups[orthogroups$Orthogroup %in% species.gain, c("Orthogroup","Eurema_hecabe")]

GenesOrtho.list<-matrix(NA, nrow = 0, ncol = 2)
GenesOrtho.list<-as.data.frame(GenesOrtho.list)
names(GenesOrtho.list)[1]="Genes"
names(GenesOrtho.list)[2]="Orthogroup"
pathways<-ortho.butterfly
names(pathways)<-c("Orthogroup","Genes")
for (i in 1:nrow(pathways)) {
  subtable<-pathways[i,]
  rcnames<-list(c(strsplit(subtable$Orthogroup[1], ',')[[1]]),c(strsplit(subtable$Genes[1], ', ')[[1]]))
  pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
  pairtable<-as.data.frame(pairtable)
  pairtable$Orthogroup<-rownames(pairtable)
  rownames(pairtable)<-1:nrow(pairtable)
  pairtable<-as.data.frame(pairtable)
  pairtable.new<-pairtable %>% gather(Genes, pair, c(1:ncol(pairtable)-1))
  pairtable.new<-pairtable.new[,c(1:2)]
  count.og<-as.numeric(base.change.tab$`Eurema_hecabe<0>`[base.change.tab$FamilyID == subtable$Orthogroup[1]])
  pairtable.new<-pairtable.new[1:count.og,]
  GenesOrtho.list<-rbind(GenesOrtho.list, pairtable.new)
}
GenesOrtho.list<-subset(GenesOrtho.list,
                        is.na(GenesOrtho.list$Genes)==FALSE,
                        select = c("Genes", "Orthogroup"))
GenesOrtho.list<-unique(GenesOrtho.list)


write.table(GenesOrtho.list,
            file = "butterfly.gain.txt",
            sep = '\t', row.names = FALSE,
            quote = FALSE)
rm(pairtable, pairtable.new, rcnames, subtable, pathways)

##function enrich
GenesOrtho.list$Groups<-"Ehe.gain"
#read in orthogroup functions
##KEGG
kegg2name<-read.delim(paste(opt$anno, "kegg2name.txt", sep = "/"), 
                      sep = "\t", colClasses = "character")
kegg2name.nh<-read.delim(paste(opt$anno, "kegg2name-nonhuman.txt", sep = "/"), 
                         sep = "\t", colClasses = "character")
kegg2ont<-read.delim(paste(opt$anno,"kegglevel.AC.txt", sep = "/"),
                     sep = "\t", colClasses = "character")
names(kegg2ont)[2]="ID"

Orthogroups.KEGG<-read.delim("GenesKEGGpair.1v1",header = TRUE, sep = " ")

##KOG
kog2name<-read.delim(paste(opt$anno, "kog2name.txt", sep = "/"), 
                     sep = "\t", colClasses = "character")
Orthogroups.KOG<-read.delim("GenesKOGpair.1v1", 
                            header = TRUE, sep = " ")
##GO
go2name<-read.delim(paste(opt$anno, "go2name.txt", sep = "/"), 
                    sep = "\t", colClasses = "character",
                    header = FALSE)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"
Orthogroups.GO<-read.delim("GenesGOpair.1v1", 
                           header = TRUE, sep = " ")


#GO
{
  GO.ehe<-compareCluster(Genes ~ Groups, 
                              data = GenesOrtho.list, 
                              fun = 'enricher',
                              TERM2GENE = Orthogroups.GO,
                              TERM2NAME = go2name,
                              pvalueCutoff = 1,
                              pAdjustMethod = "BH",
                              qvalueCutoff = 1,
                              minGSSize = 1,
                              maxGSSize = 1000000)
  GO.ehe<-as.data.frame(GO.ehe)
  
  #Plot
  plotin<-rbind(GO.ehe)
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[names(plotinsep)=="ID"]<-"goClass"
  go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)
  
  plotdata<-merge(plotinsep, go2name, by = c("goClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  
  write.table(plotdata, file = "gain/GOenrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  level3_terms<-read.delim(paste(opt$anno, "GOlevel3.txt", sep = "/"), header = TRUE)
  level4_terms<-read.delim(paste(opt$anno, "GOlevel4.txt", sep = "/"), header = TRUE)
  ##all GO
  plotdata1<-subset(plotdata, ratio2 >=0.1 & p.adjust < 0.2)
  a<-ggplot(plotdata[which(plotdata$goClass %in% plotdata1$goClass & plotdata$ratio2 > 0.1),], 
            aes(x = ratio2, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Ratio", x = "")+
    geom_point(shape = 19)+scale_size_area()+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 11))
  a
  ggsave("gain/GO.bygroup.P20.tiff", width = 12, height = 30, units = "in", dpi = 300, limitsize = FALSE)
  ggsave("gain/GO.bygroup.P20.png", width = 12, height = 30, units = "in", dpi = 300, limitsize = FALSE)
  
  ##Level 4  
  
  {

    topgroups<-plotdata1[which(plotdata1$p.adjust < 0.2),] %>% group_by(Cluster) %>% top_n(20, ratio1)
    plotdata1$ONTOLOGY<-gsub("Biological Process", "Biological\nProcess",plotdata1$ONTOLOGY)
    plotdata1$ONTOLOGY<-gsub("Molecular Function", "Molecular\nFunction",plotdata1$ONTOLOGY)
    plotdata1$ONTOLOGY<-gsub("Cellular Component", "Cellular\nComponent",plotdata1$ONTOLOGY)
    a<-ggplot(subset(plotdata1, ratio1 > 0 & p.adjust < 1 & goClass %in% topgroups$goClass), 
              aes(x = ratio2, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "", size = "Ratio")+
      geom_point(shape = 19)+scale_size_area()+
      scale_x_continuous(limits = c(0,0.25))+
      scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
      guides(size = guide_legend(order = 1))+
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
            axis.title.x = element_text(colour = "white"), 
            axis.title.y = element_text(colour = "white"), 
            axis.text = element_text(colour = "black", size = 11), 
            axis.text.x = element_text(size = 11), 
            axis.text.y = element_text(size = 11),
            legend.title = element_text(size = 11), 
            plot.title = element_text(size = 12))+
      facet_grid(ONTOLOGY~., scales = "free", space = "free")+
      theme(strip.text = element_text(size = 11))
    a
    ggsave("gain/GO.L4.top20.tiff", width = 8, height = 6, units = "in", dpi = 300, limitsize = FALSE)
    ggsave("gain/GO.L4.top20.png", width = 8, height = 6, units = "in", dpi = 300, limitsize = FALSE)
    
  }
}

#KEGG
{
  KEGG.ehe<-compareCluster(Genes ~ Groups, 
                                data = GenesOrtho.list, 
                                fun = 'enricher',
                                TERM2GENE = Orthogroups.KEGG,
                                TERM2NAME = kegg2name,
                                pvalueCutoff = 1,
                                pAdjustMethod = "BH",
                                qvalueCutoff = 1,
                                minGSSize = 1,
                                maxGSSize = 200000)
  KEGG.ehe<-as.data.frame(KEGG.ehe)
  
  ##Plot
  plotin<-rbind(KEGG.ehe)
  
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotdata<-subset(plotinsep, is.na(Description) == FALSE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  
  
  plotdata<-subset(plotdata, is.na(Description) == FALSE & p.adjust <= 1 & ID != "ko04013")
  plotdata<-merge(plotdata, kegg2ont, by = "ID", all.x = TRUE)
  unique(plotdata$Group)
  write.table(plotdata, file = "gain/KEGG.enrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  plotdata1<-plotdata[plotdata$Count>=2 & is.na(plotdata$ONTOLOGY)==F & plotdata$ONTOLOGY != "Human Diseases",]
  plotdata1$ONTOLOGY<-gsub("Environmental Information Processing", "Environmental\nInformation\nProcessing", plotdata1$ONTOLOGY)
  plotdata1$ONTOLOGY<-gsub("Genetic Information Processing", "Genetic\nInformation\nProcessing", plotdata1$ONTOLOGY)
  plotdata1$ONTOLOGY<-gsub("Organismal Systems", "Organismal\nSystems", plotdata1$ONTOLOGY)
  plotdata1$ONTOLOGY<-gsub("Cellular Processes", "Cellular\nProcesses", plotdata1$ONTOLOGY)
  
  topgroups<-plotdata1 %>% group_by(Cluster) %>% top_n(20, ratio1*Count)
  a<-ggplot(subset(plotdata1, ratio1 > 0 & ID %in% topgroups$ID ), 
            #  a<-ggplot(subset(plotdata, ratio1 >= 0 & p.adjust < 0.2), 
            aes(x = ratio1, y = Description, size = ratio1, colour = p.adjust))+
    labs(title = "", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_x_continuous(limits = c(0,0.55), breaks = seq(0,0.55,0.1))+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 11))
  a
  ggsave("gain/KEGG.top20.tiff", width = 8, height = 8, units = "in", dpi = 300, limitsize = FALSE)
  ggsave("gain/KEGG.top20.png", width = 8, height = 8, units = "in", dpi = 300, limitsize = FALSE)
  
}

#KOG
{
  KOG.ehe<-compareCluster(Genes ~ Groups, 
                               data = GenesOrtho.list, 
                               fun = 'enricher',
                               TERM2GENE = Orthogroups.KOG,
                               TERM2NAME = kog2name,
                               pvalueCutoff = 1,
                               pAdjustMethod = "BH",
                               qvalueCutoff = 1,
                               minGSSize = 1,
                               maxGSSize = 200000)
  KOG.ehe<-as.data.frame(KOG.ehe)
  
  ##Plot
  plotin<-rbind(KOG.ehe)
  #View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"
  kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
  kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
  kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
  kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)
  plotdata<-merge(plotinsep, kog2name, by = c("kogClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  
  write.table(plotdata, file = "gain/KOGenrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  plotdata$ONTOLOGY<-gsub("Information storage and processing", "Information\nstorage and\nprocessing", plotdata$ONTOLOGY)
  plotdata$ONTOLOGY<-gsub("Cellular processes and signaling", "Cellular\nprocesses and\nsignaling",plotdata$ONTOLOGY)
  a<-ggplot(plotdata, aes(x = ratio1, y = Description, size = ratio1, colour = p.adjust))+
    labs(title = "", size = "Ratio")+
    geom_point(shape = 19)+scale_size_area()+
    scale_x_continuous(limits = c(0,0.05), breaks = seq(0,0.05,0.01))+
    scale_color_gradient(low = "#FF0000", high = "#0000FF", limits = c(0,0.2))+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "white"), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(size = 11), 
          axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 11), 
          plot.title = element_text(size = 12))+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 11))
  a
  ggsave("gain/KOG.tiff", width = 8, height = 6, units = "in", dpi = 300)
  ggsave("gain/KOG.png", width = 8, height = 6, units = "in", dpi = 300)
  
}

