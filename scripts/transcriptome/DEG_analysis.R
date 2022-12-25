setwd("D:/butterfly/DEG/")
save.image("DEG.RData")
load("DEG.RData")
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(pathview)
library(clusterProfiler)
filelist<-read.delim("filelist.txt", header = TRUE)

DEG<-as.data.frame(matrix(NA, ncol = 13))
names(DEG)<-c("Genes","sampleA","sampleB","logFC","logCPM","PValue","FDR",
              "Temp.sampleA","Sex.sampleA","Tissue.sampleA",
              "Temp.sampleB","Sex.sampleB","Tissue.sampleB")

for (i in 1:nrow(filelist)) {
  DEG.sub<-read.delim(file = filelist$file[i], header = TRUE)
  DEG.sub$Genes<-row.names(DEG.sub)
  row.names(DEG.sub)<-1:nrow(DEG.sub)
  #DEG.sub<-subset(DEG.sub, abs(logFC)>1 & FDR < 0.05)
  DEG.sub$sampleA<-gsub("Ehe_24h_","", DEG.sub$sampleA)
  DEG.sub$sampleB<-gsub("Ehe_24h_","", DEG.sub$sampleB)
  DEG.sub$Temp.sampleA<-substr(DEG.sub$sampleA,1,2)
  DEG.sub$Temp.sampleA<-paste0(DEG.sub$Temp.sampleA, "°C")
  DEG.sub$Sex.sampleA<-substr(DEG.sub$sampleA,3,3)
  DEG.sub$Tissue.sampleA<-substr(DEG.sub$sampleA,5,5)
  DEG.sub$Temp.sampleB<-substr(DEG.sub$sampleB,1,2)
  DEG.sub$Temp.sampleB<-paste0(DEG.sub$Temp.sampleB, "°C")
  DEG.sub$Sex.sampleB<-substr(DEG.sub$sampleB,3,3)
  DEG.sub$Tissue.sampleB<-substr(DEG.sub$sampleB,5,5)
  DEG<-rbind(DEG, DEG.sub)
}

rm(DEG.sub)
rm(i, filelist)

DEG<-DEG[is.na(DEG$sampleA)==F,]
write.table(DEG, file = "DEG.allpairs.txt",
            sep = "\t", quote = F, row.names = T)

DEG<-read.delim("DEG.allpairs.txt", header = T)
DEG$Genes<-row.names(DEG)
DEG$change<-NA
DEG$change[DEG$logFC<0]<-"Down"
DEG$change[DEG$logFC>0]<-"Up"
DEG$Groups<-paste0(DEG$sampleA,"\n VS \n", DEG$sampleB)

DEG.summary<-as.data.frame(xtabs(~sampleA+sampleB+Temp.sampleA+Sex.sampleA+Tissue.sampleA+Temp.sampleB+Sex.sampleB+Tissue.sampleB,DEG.bk))
DEG.summary<-DEG.summary[DEG.summary$Freq>0,]
DEG.summary$sampleA1=paste0(DEG.summary$Sex.sampleA, DEG.summary$Tissue.sampleA,DEG.summary$Temp.sampleA)
DEG.summary$sampleB1=paste0(DEG.summary$Sex.sampleB, DEG.summary$Tissue.sampleB,DEG.summary$Temp.sampleB)
write.table(DEG.summary, "DEG.summary.txt", row.names = F, quote = F, sep = "\t")

##Functional annotation
addline_format <- function(x,...){ 
  gsub('_','\n',x) 
} 

GenesKOGpair.1v1<-read.delim("D://butterfly/eggnog/KOG.1v1.txt", header = TRUE)
GenesKEGGpair.1v1<-read.delim("D://butterfly/eggnog/KEGG.1v1.txt", header = TRUE)
GenesGOpair.1v1<-read.delim("D://butterfly/eggnog/GO.1v1.txt", header = TRUE)

kog2name<-read.delim("D://3enrichment/kog2name.txt", header = TRUE)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

kegg2name<-read.delim("D://3enrichment/kegg2name.txt")

kegg2ont<- read.delim("D:/3enrichment/kegglevel.AC.txt", 
                      sep = "\t", colClasses = "character")
names(kegg2ont)[2]="ID"

go2name<-read.delim("D:/3enrichment/go2name.txt", 
                    sep = "\t", colClasses = "character",
                    header = FALSE)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"
go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)

DEG.bk<-DEG
DEG.sum.sub<-DEG.summary[DEG.summary$Freq>100,]
DEG<-DEG[DEG$sampleA %in% DEG.sum.sub$sampleA & DEG$sampleB %in% DEG.sum.sub$sampleB, ]
####KOG enrich
{
  KOG.all.1<-compareCluster(Genes ~ Groups+change, 
                            data = DEG, 
                            fun = 'enricher',
                            TERM2GENE = GenesKOGpair.1v1,
                            TERM2NAME = kog2name,
                            pvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            qvalueCutoff = 1,
                            minGSSize = 1,
                            maxGSSize = 200000)
  
  ##Plot
  plotin<-as.data.frame(KOG.all.1)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[4]<-"kogClass"
  
  plotdata<-merge(plotinsep, kog2name, by = c("kogClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = "KOGenrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  list(plotdata$Group)
  a<-ggplot(plotdata, aes(x = change, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene ratio", colour = "Adjust P value", xlab = "", ylab = "")+
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
    facet_grid(ONTOLOGY~Groups, scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave("KOG.tiff", width = 12, height = 8, units = "in", dpi = 300)
  ggsave("KOG.png", width = 12, height = 8, units = "in", dpi = 300)
}


####KEGG enrich
{
  KEGG.all.1<-compareCluster(Genes ~ Groups+change, 
                             data = DEG, 
                             fun = 'enricher',
                             TERM2GENE = GenesKEGGpair.1v1,
                             TERM2NAME = kegg2name,
                             pvalueCutoff = 1,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 1,
                             minGSSize = 1,
                             maxGSSize = 200000)
  
  ##Plot
  plotin<-as.data.frame(KEGG.all.1)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  plotdata<-subset(plotinsep, is.na(Description) == FALSE & p.adjust < 1)
  
  plotdata<-merge(plotdata, kegg2ont, by = "ID", all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = "KEGG.enrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  plotdata$ONTOLOGY<-gsub(" ", "\n", plotdata$ONTOLOGY)
  list(plotdata$Group)
  
  plotdata1<-subset(plotdata, Groups != "25F_B\n VS \n25M_B" & Groups != "30F_B\n VS \n30M_B")
  a<-ggplot(subset(plotdata1, ratio2>0 & p.adjust< 0.2 & ONTOLOGY != "Human\nDiseases"),
            aes(x = change, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene count", xlab = "Gene Ratio")+
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
          plot.title = element_text(size = 12),
          strip.text.x = element_text(size = 11, angle = 0),
          strip.text.y = element_text(size = 9, angle = 0))+
    facet_grid(ONTOLOGY~Groups, scales = "free", space = "free")
    
  a
  ggsave("KEGG.temperature.H2B.tiff", width = 15, height = 20, units = "in", dpi = 300)
  ggsave("KEGG.temperature.H2B.png", width = 15, height = 20, units = "in", dpi = 300)
  
  plotdata1<-subset(plotdata, Groups == "25F_B\n VS \n25M_B" | Groups == "30F_B\n VS \n30M_B")
  a<-ggplot(subset(plotdata1, ratio2>0 & p.adjust< 0.2 & ONTOLOGY != "Human\nDiseases"),
            aes(x = change, y = Description, size = ratio2, colour = p.adjust))+
    labs(title = "", size = "Gene count", xlab = "Gene Ratio")+
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
          plot.title = element_text(size = 12),
          strip.text.x = element_text(size = 11, angle = 0),
          strip.text.y = element_text(size = 9, angle = 0))+
    facet_grid(ONTOLOGY~Groups, scales = "free", space = "free")
  
  a
  ggsave("KEGG.temperature.F2M.tiff", width = 10, height = 8, units = "in", dpi = 300)
  ggsave("KEGG.temperature.F2M.png", width = 10, height = 8, units = "in", dpi = 300)
}

####GO enrich
{
  GO.all.1<-compareCluster(Genes ~ Groups+change, 
                           data = DEG, 
                           fun = 'enricher',
                           TERM2GENE = GenesGOpair.1v1,
                           TERM2NAME = go2name,
                           pvalueCutoff = 1,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 1,
                           minGSSize = 1,
                           maxGSSize = 2000000)
  
  #Plot
  plotin<-as.data.frame(GO.all.1)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[4]<-"goClass"
  go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
  go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)
  
  plotdata<-merge(plotinsep, go2name, by = c("goClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  #test<-plotdata[which(plotdata$goClass %in% GOl5$goClass),]
  
  write.table(plotdata, file = "GOenrich.txt",
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  GO.L4<-read.delim("D:/3enrichment/GOlevel4.txt")
  plotdata<-plotdata[-which(plotdata$goClass %in% GO.L4$goClass),]
  
  GO.L5<-read.delim("D:/3enrichment/GOlevel5.txt")
  plotdata<-plotdata[-which(plotdata$goClass %in% GO.L5$goClass),]
  
  GO.L6<-read.delim("D:/3enrichment/GOlevel6.txt")
  plotdata<-plotdata[-which(plotdata$goClass %in% GO.L6$goClass),]
  
  GO.L7<-read.delim("D:/3enrichment/GOlevel7.txt")
  plotdata<-plotdata[which(plotdata$goClass %in% GO.L7$goClass),]
  
  ##all GO
  plotdata1<-subset(plotdata, Groups != "25F_B\n VS \n25M_B" & Groups != "30F_B\n VS \n30M_B")
    a<-ggplot(subset(plotdata1, ratio2>0.05 & p.adjust<0.2),
              aes(x = change, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "", size = "Gene count", xlab = "Gene Ratio")+
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
      facet_grid(ONTOLOGY~Groups, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 8))
    a
    
    ggsave("GO.temperature.H2B.tiff", width = 10, height = 12, units = "in", dpi = 300)
    ggsave("GO.temperature.H2B.png", width = 10, height = 12, units = "in", dpi = 300)
    
    plotdata1<-subset(plotdata, Groups == "25F_B\n VS \n25M_B" | Groups == "30F_B\n VS \n30M_B")
    a<-ggplot(subset(plotdata1, ratio2>0.01 & p.adjust<0.2 & ONTOLOGY != "Cellular Component"),
              aes(x = change, y = Description, size = ratio2, colour = p.adjust))+
      labs(title = "", size = "Gene count", xlab = "Gene Ratio")+
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
      facet_grid(ONTOLOGY~Groups, scales = "free", space = "free")+
      theme(strip.text = element_text(size = 8))
    a
    ggsave("GO.temperature.F2M.tiff", width = 8, height = 8, units = "in", dpi = 300)
    ggsave("GO.temperature.F2M.png", width = 8, height = 8, units = "in", dpi = 300)
}

##pathview mapping
DEG.summary<- read.delim("DEG.summary.txt", header = TRUE)
DEG.sum.sub<-subset(DEG.summary, pairtype == "temp")
DEG<-DEG.bk
DEG.bk$Genes<-row.names(DEG.bk)
DEG.sum.sub<-DEG.sum.sub[,1:2]
DEG.sex<-merge(DEG.sum.sub, DEG.bk, by = c("sampleA", "sampleB"), all.x = TRUE)
DEG.sex$Genes<-DEG.sex(DEG)
GenesKOpair.1v1<-read.delim("D:/butterfly/eggnog/KO.1v1.txt", header = TRUE)
DEG.sex<-merge(DEG.sex, GenesKOpair.1v1, all.x = TRUE, by = "Genes")
DEG.sex<-DEG.sex[is.na(DEG.sex$ko)==F,]
KEGG2koID<-read.delim("D:/3enrichment/kegg2koID.txt",header = TRUE)
DEG.sex<-merge(DEG.sex, KEGG2koID, by = "ko", all.x = TRUE)
pathways<-as.data.frame(unique(DEG.sex$KEGG))
names(pathways)[1]="KEGG"
pathways<-merge(pathways, kegg2name, by = "KEGG", all.x = TRUE)
names(pathways)[1]="pathways"
write.table(pathways, "pathways.txt", row.names = F, quote = F, sep = "\t")

DEG.matrix<-as.data.frame(unique(DEG$Genes))
names(DEG.matrix)[1]="Genes"

for (i in 1:nrow(DEG.sum.sub)) { 
  DEG.sub<-subset(DEG, sampleA == DEG.sum.sub$sampleA[i] & sampleB == DEG.sum.sub$sampleB[i],
                 select = c("Genes", "logFC"))
  names(DEG.sub)[2]=paste0(DEG.sum.sub$sampleA[i],".", DEG.sum.sub$sampleB[i])
  DEG.sub<-unique(DEG.sub)
  DEG.matrix<-merge(DEG.matrix, DEG.sub, all.x = TRUE, by = "Genes")
}

DEG.matrix<-merge(DEG.matrix, GenesKOpair.1v1, all.x = TRUE, by = "Genes")
DEG.matrix<-DEG.matrix[,-1]
DEG.ko.matrix<-as.data.frame(unique(DEG.matrix$ko))
names(DEG.ko.matrix)[1]="ko"

for (i in 1:(ncol(DEG.matrix)-2)){
  DEG.matrix.sub<-DEG.matrix[,c(i,9)]
  names(DEG.matrix.sub)[1]="FC"
  DEG.ko.matrix.sub<-DEG.matrix.sub %>% group_by(ko) %>% summarise(sum = sum(FC))
  names(DEG.ko.matrix.sub)[2]=names(DEG.matrix)[i]
  DEG.ko.matrix<-merge(DEG.ko.matrix, DEG.ko.matrix.sub, all.x = TRUE, by = "ko")
}
rm(DEG.matrix.sub, DEG.ko.matrix.sub)

DEG.ko.matrix<-DEG.ko.matrix[is.na(DEG.ko.matrix$ko)==F,]
row.names(DEG.ko.matrix)<-DEG.ko.matrix$ko
DEG.ko.matrix<-DEG.ko.matrix[,-1]
DEG.ko.matrix<-DEG.ko.matrix[,c(1,2,3,6,4,7,5,8)]
write.table(DEG.ko.matrix, "DEG.ko.matrix.sex.txt",
            sep = "\t", quote = F, row.names = T)
setwd("pv/")
for (i in 1:nrow(pathways)) {
  pathwayid<-pathways$pathways[i]
  pathwayid<-"ko04150"
  pathview(gene.data = DEG.ko.matrix,
           pathway.id = pathwayid, 
           species = "ko", 
           gene.idtype = "KEGG", 
           limit = list(gene = 1), 
           bins = list(gene=10), 
           multi.state = TRUE, 
           na.col="transparent", 
           out.suffix = "sex_bias")
}
