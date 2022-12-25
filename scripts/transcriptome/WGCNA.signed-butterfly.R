#sink("WGCNA-signed.cold.txt", type=c("output", "message"))
sink("WGCNA-signed.head.txt", type=c("output", "message"))
setwd("D:/butterfly/WGCNA/")
##WGCNA
library(WGCNA)
library(psych)
library(tidyr)
library(Hmisc)
library(dplyr)
library(scatterplot3d)

####read in expression
exp.table<-read.delim(file = "RSEM.gene.TMM.EXPR.matrix", 
                      header = TRUE)
row.names(exp.table)=exp.table$X
exp.mRNA.sum<-exp.table[,c(1,2)]
exp.table<-exp.table[,2:ncol(exp.table)]

#log2 transformation
exp.table<-log2(exp.table+1)

####read in trait
trait<-read.delim(file = "trait.txt", header = TRUE)
trait.bk1<-trait
rownames(trait) = trait$ID
trait <- subset(trait, select = -c(ID))
trait.bk<-trait
trait.heat <- trait[which(is.na(trait$Temperature.Heat) == FALSE),-3]
trait.cold <- trait[which(is.na(trait$Temperature.Cold) == FALSE),-4]
trait.head <- trait[which(trait$Head ==1),c(1,3,4)]
trait.body <- trait[which(trait$Body ==1),c(1,3,4)]

##Treatment
#trait<-trait.cold
#treatment="cold"

treatment="all"
Rdata="WGCNA.all.RData"

trait<-trait.head
treatment="head"
Rdata="WGCNA.head.RData"

trait<-trait.body
treatment="body"
Rdata="WGCNA.body.RData"

exp.mRNA<-exp.table[, which(names(exp.table) %in% row.names(trait))]

##PCA
pca.matrix<-prcomp(t(exp.mRNA))
pca.eigenvalues <-pca.matrix$sdev^2
####% of explanation
pca.eigenvalues <- tibble(PC = factor(1:length(pca.eigenvalues)), 
                         variance = pca.eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))

p<-pca.eigenvalues %>% 
  ggplot(aes(x = PC))+
  geom_col(aes(y = pct))+
  geom_line(aes(y = pct_cum, group = 1))+ 
  geom_point(aes(y = pct_cum))+
  scale_y_continuous(breaks = seq(0,100,10))+
  labs(x = "Principal component", y = "Fraction variance explained")
p  
ggsave(paste0(treatment,".PCA.variance.sum.png"), 
       width = 6, height = 4.5, units = "in", dpi = 300)

####distribution
pca.scores <- pca.matrix$x
pca.scores <- pca.scores %>% 
  as_tibble(rownames = "sample")
names(trait.bk1)[1]="sample"
trait.bk1$Part[which(is.na(trait.bk1$Head)==FALSE & trait.bk1$Sexuality == -1)]<-"Female head"
trait.bk1$Part[which(is.na(trait.bk1$Head)==FALSE & trait.bk1$Sexuality == 1)]<-"Male head"
trait.bk1$Part[which(is.na(trait.bk1$Body)==FALSE & trait.bk1$Sexuality == -1)]<-"Female body"
trait.bk1$Part[which(is.na(trait.bk1$Body)==FALSE & trait.bk1$Sexuality == 1)]<-"Male body"

pca.scores<-merge(pca.scores, trait.bk1, by = "sample", all.x = TRUE)
####distribution
pca.scores <- pca.matrix$x
pca.scores <- pca.scores %>% 
  as_tibble(rownames = "sample")
names(trait.bk1)[1]="sample"

trait.bk1$Sex<-gsub("-1","Female", trait.bk1$Sex)
trait.bk1$Sex<-gsub("1","Male", trait.bk1$Sex)

pca.scores<-merge(pca.scores, trait.bk1, by = "sample", all.x = TRUE)
pca.scores$Temperature<-paste0(pca.scores$Temperature, "°C")
p<-pca.scores %>% 
  ggplot(aes(x = PC1, y = PC2, color = factor(Part))) +
  geom_point(aes(#fill = factor(Temperature),
    shape = factor(Temperature)),
    alpha = 0.6, size = 3)+
  #geom_text(aes(label = sample))+
  #scale_shape_manual(values = c(19,17), limits = c("Male","Female"))+
  labs(x = paste("Principal Component 1 (explained variance = ", round(pca.eigenvalues$pct[1],2), "%)", sep = ""),
       y = paste("Principal Component 2 (explained variance = ", round(pca.eigenvalues$pct[2],2), "%)", sep = ""),
       shape = "Temperature", color = "Part")
p
ggsave(paste0(treatment,".PCA.PC1-PC2.png"), 
       width = 6, height = 4.5, units = "in", dpi = 300)

p<-pca.scores %>% 
  ggplot(aes(x = PC2, y = PC3, color = factor(Part))) +
  geom_point(aes(#fill = factor(Temperature),
    shape = factor(Temperature)),
    alpha = 0.6, size = 3)+
  #geom_text(aes(label = sample))+
  #scale_shape_manual(values = c(19,17), limits = c("Male","Female"))+
  labs(x = paste("Principal Component 2 (explained variance = ", 
                 round(pca.eigenvalues$pct[2],2), "%)", sep = ""),
       y = paste("Principal Component 3 (explained variance = ", 
                 round(pca.eigenvalues$pct[3],2), "%)", sep = ""),
       shape = "Temperature", color = "Part")
p
ggsave(paste0(treatment,".PCA.PC2-PC3.png"),
       width = 6, height = 4.5, units = "in", dpi = 300)

p<-pca.scores %>% 
  ggplot(aes(x = PC1, y = PC3, color = factor(Part))) +
  geom_point(aes(#fill = factor(Temperature),
    shape = factor(Temperature)),
    alpha = 0.6, size = 3)+
  #geom_text(aes(label = sample))+
  #scale_shape_manual(values = c(19,17), limits = c("Male","Female"))+
  labs(x = paste("Principal Component 1 (explained variance = ", 
                 round(pca.eigenvalues$pct[1],2), "%)", sep = ""),
       y = paste("Principal Component 3 (explained variance = ", 
                 round(pca.eigenvalues$pct[3],2), "%)", sep = ""),
       shape = "Temperature", color = "Part")
p
ggsave(paste0(treatment,".PCA.PC1-PC3.png"), 
              width = 6, height = 4.5, units = "in", dpi = 300)

##sample correlation
d<-cor(exp.mRNA, method = "pearson", use = "na.or.complete")
h<-hclust(as.dist(1-d), method = "complete")

png(filename = paste(treatment,".sample_pearson_tree.png", sep = ""), 
    res = 300, units = "in", width = 8, height = 6)
plot(varclus(d, similarity="pearson", type = "similarity.matrix"), hang = -1)
dev.off()

##filter by MAD > 0
{
  exp.mRNA.sum$max = apply(exp.mRNA, 1, max)
  exp.mRNA.sum$mean = apply(exp.mRNA, 1, mean)
  exp.mRNA.sum<-exp.mRNA.sum[,c(3:4)]
  exp.mRNA.sum$mad = apply(exp.mRNA, 1, mad)
  exp.mRNA.sum$TMMcut = rowSums(as.matrix(exp.mRNA) > log2(10))
  exp.mRNA.sumsub<-subset(exp.mRNA.sum, 
                          mad > 0 & TMMcut >= 3)
}

##filter matrix
exp.mRNA<-exp.mRNA[which(rownames(exp.mRNA) %in% rownames(exp.mRNA.sumsub)),]
rm(exp.mRNA.sum, exp.mRNA.sumsub)

##Sample correlation
d<-cor(exp.mRNA, method = "pearson", use = "na.or.complete")
h<-hclust(as.dist(1-d), method = "complete")

png(filename = paste(treatment, ".sample_pearson_tree.filter.png", sep = ""),
    res = 300, units = "in", width = 8, height = 6)
plot(varclus(d, similarity="pearson", type = "similarity.matrix"), hang = -1)
dev.off()
rm(d,h)

##enter file
write.table(exp.mRNA, 
            file = paste(treatment, ".filter.matrix.txt", sep = ""),
            quote = F, row.names = T, sep = "\t")
exp.mRNA<-read.delim(paste(treatment, ".filter.matrix.txt", sep = ""),
                     header = T)
#transform matrix
WGCNA.data <- t(exp.mRNA)
nSamples <- nrow(WGCNA.data)

allowWGCNAThreads(6)
save.image(Rdata)

###########################################################
####signed
###########################################################

networkType="signed"
cortype="bicor"

####pickSoftThreshold
powers = c(seq(from = 1, to=30, by=1))
sft = pickSoftThreshold(WGCNA.data, 
                        powerVector = powers, 
                        networkType = networkType, 
                        verbose = 5)

png(file = paste(treatment, ".", cortype, ".SoftThreshold.png", sep = ""),
    width = 6, height = 6*0.75, units = "in", res = 300)
SoftThreshold<-plot(sft$fitIndices[,1], 
                    -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
                    xlab='Soft Threshold (power)', 
                    ylab='Scale Free Topology Model Fit, signed R^2', 
                    type = 'n', main = paste('Scale independence'))+
  text(sft$fitIndices[,1],
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=1,col='red')+
  abline(h=0.85,col='red')
dev.off()

####pick the first power that reach 0.85
beta=10
cat(paste("beta=",beta, "\n", sep = ""))  

save.image(Rdata)

####Mean Connectivity
png(file = paste(treatment, ".", cortype, ".meanconnectivity.png", sep = ""),
    width = 6, height = 6*0.75, units = "in", res = 300)
meanconnectivity<-plot(sft$fitIndices[,1], sft$fitIndices[,5], 
                       xlab="Soft Threshold (power)", ylab="Mean Connectivity", 
                       type = 'n', main = paste("Mean connectivity:", cortype, sep=""))+         
  text(sft$fitIndices[,1], sft$fitIndices[,5], 
       labels=powers, cex=1, col="red")+
  abline(h=100,col='red')
dev.off()
save.image(Rdata)

####network construction
net=blockwiseModules(WGCNA.data, power = beta,
                     maxBlockSize = 30000, 
                     minModuleSize = 30,
                     #mergeCutHeight = 0.20, #head
                     mergeCutHeight = 0.6, #body
                     corType = cortype, 
                     networkType = networkType,
                     TOMType = networkType, 
                     saveTOMs = TRUE, 
                     nThreads = 6,
                     maxPOutliers = 0.05,
                     verbose = 5)

save.image(Rdata)

table(net$colors)
write.table(table(net$colors), 
            file = paste(treatment, ".", cortype, ".net.colors.txt", sep = ""),
            sep = '\t', row.names = FALSE,
            quote = FALSE)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
png(file = paste(treatment, ".", cortype, ".dendrograms.png", sep = ""),
    width = 6, height = 6, units = "in", res = 300)
dendrograms<-plotDendroAndColors(net$dendrograms[[1]], 
                                 moduleColors[net$blockGenes[[1]]], 
                                 "Module colors", dendroLabels = FALSE, hang = 0.03, 
                                 addGuide = TRUE, guideHang = 0.05)         
dev.off()
save.image(Rdata)

####correlation among modules
MEs = net$MEs
MEs_col = MEs
MEs_col = orderMEs(MEs_col)
png(file = paste(treatment, ".", cortype, ".Eigengene.png", sep = ""),
    width = 20, height = 20*2, units = "in", res = 300)
Eigengene<-plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                                 marDendro = c(3,3,2,4), marHeatmap = c(3,4,2,2), 
                                 plotDendrograms = T, xLabelsAngle = 90,
                                 colorLabels = T, signed = T,
                                 plotAdjacency = T,
                                 printAdjacency = T)
dev.off()
save.image(Rdata)

####Trait
#trait<-trait.body
#trait$Head[is.na(trait$Head)==T]<-0
#trait$Body[is.na(trait$Body)==T]<-0
#trait<-trait[,-2]
if (cortype == "bicor") {
  moduleTraitCor.bicor = bicorAndPvalue(MEs_col, trait)
  moduleTraitCor = moduleTraitCor.bicor$bicor
  moduleTraitPvalue = moduleTraitCor.bicor$p
  ind.signif=which(moduleTraitCor.bicor$p<2 & abs(moduleTraitCor.bicor$bicor) >= 0,
                   arr.ind = T)
  bicor.signif<-moduleTraitCor.bicor$bicor[ind.signif]
  p.signif<-moduleTraitCor.bicor$p[ind.signif]
  bicor.signif.module<-dimnames(moduleTraitCor.bicor$bicor)[[1]][ind.signif[,1]]
  bicor.signif.trait<-dimnames(moduleTraitCor.bicor$bicor)[[2]][ind.signif[,2]]
  bicor.signif.filter<-cbind(bicor.signif.module, bicor.signif.trait,bicor.signif,p.signif)
} else {
  moduleTraitCor = cor(MEs_col, trait, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(MEs_col))
}

textMatrix = paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep = ""); 

png(file = paste(treatment, ".", cortype, ".Module-trait.relationships.png",sep = ""),
    width = 8, height = 6, units = "in", res = 300)
ModuleTrait.relationheatmap<-
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(trait),
                 xSymbols = names(trait),
                 yLabels = names(MEs_col),
                 ySymbols = names(MEs_col),
                 colorLabels = NULL,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = TRUE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 xLabelsAngle = 0,
                 #xColorWidth = 3, 
                 xLabelsAdj = 0.5,
                 main = paste("Module-trait relationships"))
dev.off()

save.image(Rdata)

##TOM plot
load(Rdata)
if (length(net$TOMFiles) == 1) {
  load(net$TOMFiles[1], verbose=T)
} else {
  TOM = TOMsimilarityFromExpr(WGCNA.data, 
                              power = beta,
                              corType = cortype, 
                              networkType = networkType,
                              TOMType = networkType,
                              maxPOutliers = 0.05,
                              nThreads = 6)
}

TOM <- as.matrix(TOM)
#dissTOM = 1-TOM
#plotTOM = dissTOM^beta

#png(file = paste(treatment,".", cortype,".TOMplot.large.png", sep = ""),
#    width = 20, height = 20, units = "in", res = 300)
#TOMplot<-TOMplot(plotTOM, net$dendrograms,
#                 moduleColors,  
#                 main = "Network heatmap plot, all genes")
#dev.off()

####Export to cytoscape
probes = colnames(WGCNA.data)

##by color
modulename<-as.data.frame(names(MEs))
modulename$`names(MEs)`<-gsub("ME", "", modulename$`names(MEs)`)
dimnames(TOM)<-list(probes, probes)
for (i in 1:nrow(modulename)) {
  modules= modulename[i,1]
  inmodules=is.finite(match(moduleLabels, modules))
  modProbes=probes[inmodules]
  modTOM=TOM[modProbes,modProbes]
  cyt = exportNetworkToCytoscape(modTOM, 
                                 edgeFile = paste(treatment,".", cortype, ".", modules, ".edges.txt", sep=""), 
                                 nodeFile = paste(treatment,".", cortype, ".", modules, ".nodes.txt", sep=""), 
                                 weighted = TRUE, threshold = 0.2, nodeNames = modProbes, 
                                 nodeAttr = moduleLabels[inmodules])
}

#cyt = exportNetworkToCytoscape(TOM, 
                               edgeFile = paste(treatment,".", cortype, ".all", ".edges.txt", sep=""), 
                               nodeFile = paste(treatment,".", cortype, ".all", ".nodes.txt", sep=""), 
                               weighted = TRUE, threshold = 0.5, nodeNames = probes, 
                               nodeAttr = moduleColors)

##Gene color table
gene.color<-cbind(probes, moduleLabels)
gene.color<-as.data.frame(gene.color)

##Interested module to phenotype
##High correlation gene in the module
if (cortype == "bicor"){
  geneModuleMembership.bicor = bicorAndPvalue(WGCNA.data, MEs_col)
  geneModuleMembership = geneModuleMembership.bicor$bicor
  MMPvalue = geneModuleMembership.bicor$p
  geneTraitMembership.bicor = bicorAndPvalue(WGCNA.data, trait)
  geneTraitCor = as.data.frame(geneTraitMembership.bicor$bicor)
  geneTraitP = as.data.frame(geneTraitMembership.bicor$p)
} else {
  geneModuleMembership = as.data.frame(cor(WGCNA.data, MEs_col, 
                                           method = "pearson", use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  geneTraitCor = as.data.frame(cor(WGCNA.data, trait, 
                                   method = "pearson", use = "all.obs"))
  geneTraitP = as.data.frame(corPvalueStudent(as.matrix(geneTraitCor), nSamples))
}

bicor.signif.filter<-as.data.frame(bicor.signif.filter)
geneModuleMembership.bicor<-as.data.frame(geneModuleMembership.bicor)
geneTraitMembership.bicor<-as.data.frame(geneTraitMembership.bicor)

for (i in 1:nrow(bicor.signif.filter)){
  module = bicor.signif.filter$bicor.signif.module[i]
  pheno = bicor.signif.filter$bicor.signif.trait[i]
  
  #Get correlation value
  gene.module.bicor<-geneModuleMembership.bicor[, which(names(geneModuleMembership.bicor)==paste("bicor.", module, sep = ""))]
  gene.module.pvalue<-geneModuleMembership.bicor[, which(names(geneModuleMembership.bicor)==paste("p.", module, sep = ""))]
  gene.trait.bicor<-geneTraitMembership.bicor[,which(names(geneTraitMembership.bicor)==paste("bicor.", pheno, sep = ""))]
  gene.trait.pvalue<-geneTraitMembership.bicor[,which(names(geneTraitMembership.bicor)==paste("p.", pheno, sep = ""))]
  
  module.pheno.table<-cbind(gene.module.bicor, 
                            gene.module.pvalue,
                            gene.trait.bicor,
                            gene.trait.pvalue)
  module.pheno.table<-as.data.frame(module.pheno.table)
  row.names(module.pheno.table)<-row.names(geneTraitMembership.bicor)
  gene.inmodule<-gene.color[which(gene.color$moduleLabels %in% gsub("ME", "", module)), 1]
  module.pheno.table<-module.pheno.table[which(row.names(module.pheno.table) %in% gene.inmodule),]
  
  write.table(module.pheno.table,
              file = paste(treatment,".", cortype, 
                           ".Module_membership_gene_significance.", 
                           module, ".", pheno, ".bicor.txt", sep=""), 
              sep = '\t', row.names = TRUE,
              quote = FALSE)
  
  ####Trait significantly correlated genes
  x = as.numeric(as.character(module.pheno.table$gene.module.bicor))
  y = as.numeric(as.character(module.pheno.table$gene.trait.bicor))
  corExpr = parse(text = paste("cor", "(x, y ", prepComma("use = 'p'"),")"))
  cor = signif(eval(corExpr), 2)
  corp = signif(corPvalueStudent(cor, sum(is.finite(x) & is.finite(y))), 
                2)
  bicor.signif.filter$gene.cor[i]=cor
  bicor.signif.filter$gene.pvalue[i]=corp
  
  png(file = paste(treatment, ".", cortype, 
                   ".Module_membership_gene_significance.", 
                   module, ".", pheno, ".png",
                   sep = ""), 
      width = 12, height = 12, units = "in", res = 300)
  verboseScatterplot(module.pheno.table$gene.module.bicor,
                     module.pheno.table$gene.trait.bicor,
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for", pheno),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, 
                     #bg = module.pheno.table[1:10, gene.trait.pvalue]
                     )
  dev.off()
  
  moduleGenes<-row.names(module.pheno.table)
  modTOM=TOM[moduleGenes,moduleGenes]
  cyt = exportNetworkToCytoscape(modTOM, 
                                 edgeFile = paste(treatment,".", cortype, 
                                                  ".Module_membership_gene_significance.", 
                                                  module, ".", pheno, ".edges.txt", sep=""), 
                                 nodeFile = paste(treatment,".", cortype,
                                                  ".Module_membership_gene_significance.", 
                                                  module, ".", pheno, ".nodes.txt", sep=""), 
                                 weighted = TRUE, threshold = 0.2, nodeNames = modProbes, 
                                 nodeAttr = moduleLabels[inmodules])
  
}

write.table(bicor.signif.filter, 
            file = paste(treatment,".", cortype, ".bicor.signif.filter.txt", sep = ""),
            sep = '\t', row.names = FALSE,
            quote = FALSE)
rm(cor, corp, x, y, corExpr)
rm(TOM)
save.image(Rdata)

##Annotation
enrich.input<-gene.color
enrich.input<-as.data.frame(enrich.input)
write.table(enrich.input, 
            file = paste(treatment,".", cortype, ".gene.color.txt", sep = ""),
            sep = '\t', row.names = FALSE,
            quote = FALSE)
sink()
