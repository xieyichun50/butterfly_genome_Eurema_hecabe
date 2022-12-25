##WGCNA
library(WGCNA)
library(psych)
library(tidyr)
library(Hmisc)
library(dplyr)

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
rownames(trait) = trait$ID
trait <- subset(trait, select = -c(ID))
trait.heat <- trait[which(is.na(trait$Temperature.Heat) == FALSE),]
trait.cold <- trait[which(is.na(trait$Temperature.Cold) == FALSE),]

##Treatment
treatment="cold"
exp.mRNA<-exp.table[, which(names(exp.table) %in% row.names(trait.cold))]

##Treatment
#treatment="heat"
#exp.mRNA<-exp.table[, which(names(exp.table) %in% row.names(trait.heat))]

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
  exp.mRNA.sum$TMMcut = rowSums(as.matrix(exp.mRNA) > 10)
  exp.mRNA.sumsub<-subset(exp.mRNA.sum, 
                          mad > 0 & max > 1)
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
treatment="cold"
write.table(exp.mRNA, 
            file = paste(treatment, ".filter.matrix.txt", sep = ""),
            quote = F, row.names = F, sep = "\t")
exp.mRNA<-read.delim(paste(treatment, ".filter.matrix.txt", sep = ""),
                     header = T)
#transform matrix
WGCNA.data <- t(exp.mRNA)
nSamples <- nrow(WGCNA.data)

allowWGCNAThreads(6)
save.image("WGCNA.RData")

###########################################################
####Unsigned
###########################################################
#sink("WGCNA-unsigned.txt", type=c("output", "message"))
####pickSoftThreshold
#similarity.unsigned<-abs(bicor(WGCNA.data))
powers.unsigned = c(seq(from = 1, to=30, by=1))
sft.unsigned = pickSoftThreshold(WGCNA.data, 
                                 powerVector = powers.unsigned, 
                                 networkType = "unsigned" , 
                                 verbose = 5)

png(file = paste(treatment, ".SoftThreshold.unsigned.png", sep = ""),
    width = 6, height = 6*0.75, units = "in", res = 300)
SoftThreshold.unsigned<-plot(sft.unsigned$fitIndices[,1], 
                             -sign(sft.unsigned$fitIndices[,3])*sft.unsigned$fitIndices[,2], 
                             xlab='Soft Threshold (power)', 
                             ylab='Scale Free Topology Model Fit, unsigned R^2', 
                             type = 'n', main = paste('Scale independence'))+
  text(sft.unsigned$fitIndices[,1],
       -sign(sft.unsigned$fitIndices[,3])*sft.unsigned$fitIndices[,2],
       labels=powers.unsigned,cex=1,col='red')+
  abline(h=0.80,col='red')
dev.off()

####pick the first power that reach 0.80

beta.unsigned=3

cat(paste("beta.unsigned=",beta.unsigned, sep = ""))  

save.image("WGCNA.RData")

####Mean Connectivity
png(file = paste(treatment, ".meanconnectivity.unsigned.png", sep = ""),
    width = 6, height = 6*0.75, units = "in", res = 300)
meanconnectivity.unsigned<-plot(sft.unsigned$fitIndices[,1], sft.unsigned$fitIndices[,5], 
                                xlab="Soft Threshold (power)", ylab="Mean Connectivity", 
                                type = 'n', main = paste("Mean connectivity, unsigned"))+
  text(sft.unsigned$fitIndices[,1], sft.unsigned$fitIndices[,5], 
       labels=powers.unsigned, cex=1, col="red")
dev.off()
save.image("WGCNA.RData")

####network construction
net.unsigned=blockwiseModules(WGCNA.data, power = beta.unsigned,
                              maxBlockSize = 30000, 
                              minModuleSize = 30,
                              mergeCutHeight = 0.15,
                              corType = "pearson", 
                              TOMType = "unsigned", saveTOMs = TRUE, verbose = 5)

table(net.unsigned$colors)
write.table(table(net.unsigned$colors), 
            file = paste(treatment,".net.unsigned.colors.txt", sep = ""),
            sep = '\t', row.names = FALSE,
            quote = FALSE)

moduleLabels.unsigned = net.unsigned$colors
moduleColors.unsigned = labels2colors(moduleLabels.unsigned)
png(file = paste(treatment, ".dendrograms.unsigned.png", sep = ""),
    width = 6, height = 6, units = "in", res = 300)
dendrograms.unsigned<-plotDendroAndColors(net.unsigned$dendrograms[[1]], 
                                          moduleColors.unsigned[net.unsigned$blockGenes[[1]]], 
                                          "Module colors", dendroLabels = FALSE, hang = 0.03, 
                                          addGuide = TRUE, guideHang = 0.05)
dev.off()
save.image("WGCNA.RData")

####correlation among modules
MEs.unsigned = net.unsigned$MEs
MEs_col.unsigned = MEs.unsigned
MEs_col.unsigned = orderMEs(MEs_col.unsigned)
png(file = "Eigengene.unsigned.png", width = 8, height = 8*2, units = "in", res = 300)
Eigengene.unsigned<-plotEigengeneNetworks(MEs_col.unsigned, "Eigengene adjacency heatmap", 
                                          marDendro = c(3,3,2,4), marHeatmap = c(3,4,2,2), 
                                          plotDendrograms = T, xLabelsAngle = 90,
                                          colorLabels = T, signed = F,
                                          plotAdjacency = T,
                                          printAdjacency = T)
dev.off()
save.image("WGCNA.RData")

####Trait
moduleTraitCor.unsigned = cor(MEs_col.unsigned, trait, use = "p")
moduleTraitPvalue.unsigned = corPvalueStudent(moduleTraitCor.unsigned, nrow(MEs_col.unsigned))
textMatrix = paste(signif(moduleTraitCor.unsigned, 2), "\n(", signif(moduleTraitPvalue.unsigned, 1), ")", sep = ""); 

png(file = "Module-trait.relationships.png", width = 9, height = 12, units = "in", res = 300)
ModuleTrait.relationheatmap<-
  labeledHeatmap(Matrix = moduleTraitCor.unsigned,
                 xLabels = names(trait),
                 xSymbols = names(trait),
                 yLabels = names(MEs_col.unsigned),
                 ySymbols = names(MEs_col.unsigned),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = TRUE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
dev.off()

save.image("WGCNA.RData")

##High correlation gene in the module
geneModuleMembership = as.data.frame(cor(WGCNA.data, MEs_col.unsigned, 
                                         method = "pearson", use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(
  as.matrix(geneModuleMembership), nSamples))
geneTraitCor = as.data.frame(cor(WGCNA.data, trait, 
                                 method = "pearson", use = "all.obs"))
geneTraitP = as.data.frame(corPvalueStudent(
  as.matrix(geneTraitCor), nSamples))
module = "blue"
pheno = "Temperature.Heat"
modNames = substring(colnames(MEs_col.unsigned), start = 3)

module_column = match(module, modNames)
pheno_column = match(pheno,colnames(trait))

####GET GENES
inmodule=is.finite(match(moduleLabels.unsigned, module))
moduleGenes=moduleLabels.unsigned[inmodule]

####Trait significantly correlated genes
png(file = paste("Module_membership_gene_significance.", module, pheno, ".png", sep = ""), 
    width = 12, height = 12, units = "in", res = 300)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

##TOM plot
#load("WGCNA.RData")
#load(net.unsigned$TOMFiles[1], verbose=T)
#TOM.unsigned <- as.matrix(TOM)

##TOM Plot
TOM.unsigned = TOMsimilarityFromExpr(WGCNA.data, power=beta.unsigned,
                                     corType="pearson", networkType="unsigned")
TOM.unsigned <- as.matrix(TOM.unsigned)
dissTOM.unsigned = 1-TOM.unsigned
plotTOM.unsigned = dissTOM.unsigned^beta.unsigned

png(file = paste(treatment,".TOMplot.unsigned.large.png", sep = ""),
    width = 20, height = 20, units = "in", res = 300)
TOMplot.unsigned<-TOMplot(plotTOM.unsigned, net.unsigned$dendrograms,
                          moduleColors.unsigned,  
                          main = "Network heatmap plot, all genes")
dev.off()

####Export to cytoscape
probes = colnames(WGCNA.data)

##by color
modulename<-as.data.frame(names(MEs.unsigned))
modulename$`names(MEs.unsigned)`<-gsub("ME", "", modulename$`names(MEs.unsigned)`)
dimnames(TOM.unsigned)<-list(probes, probes)
for (i in 1:nrow(modulename)) {
  modules= modulename[i,1]
  inmodules=is.finite(match(moduleLabels.unsigned, modules))
  modProbes=probes[inmodules]
  modTOM=TOM.unsigned[modProbes,modProbes]
  cyt = exportNetworkToCytoscape(modTOM, edgeFile = paste(modules, ".unsigned", ".edges.txt", sep=""), 
                                 nodeFile = paste(modules, ".unsigned", ".nodes.txt", sep=""), 
                                 weighted = TRUE, threshold = 0.2, nodeNames = modProbes, 
                                 nodeAttr = moduleLabels.unsigned[inmodules])
}

cyt = exportNetworkToCytoscape(TOM.unsigned, edgeFile = paste("all", ".unsigned", ".edges.txt", sep=""), 
                               nodeFile = paste("all", ".unsigned", ".nodes.txt", sep=""), 
                               weighted = TRUE, threshold = 0.5, nodeNames = probes, 
                               nodeAttr = moduleColors.unsigned)

##Annotation
enrich.input<-cbind(probes, moduleLabels.unsigned)
enrich.input<-as.data.frame(enrich.input)
write.table(enrich.input, file = "gene.color.txt", 
            sep = '\t', row.names = FALSE,
            quote = FALSE)
