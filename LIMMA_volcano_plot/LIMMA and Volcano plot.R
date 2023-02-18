library(GEOquery)
library(limma)

gset <- getGEO("GSE110224", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))

gsms <- "0101010101010101010101010101010101"
sml <- strsplit(gsms, split="")[[1]]

ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

gs <- factor(sml)
groups <- make.names(c("Normal","Cancer"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

##################################################################################################
#                                        Volcano Plot                                            #
##################################################################################################

library(EnhancedVolcano)
res1 <- read.csv("GSE110224.top.table.csv", header=TRUE)

vol1=EnhancedVolcano(res1,
                     lab = res1$Gene.symbol,
                     x = 'logFC',
                     y = 'adj.P.Val',
                     selectLab = c('CXCL8','CEMIP','MMP7','CA4','ADH1C','GUCA2A',
                                   'GUCA2B','ZG16','CLCA4','MS4A12','CLDN1'),
                     xlab = bquote(~Log[2]~ 'fold change'),
                     pCutoff = .05,
                     FCcutoff = 1.0,
                     pointSize = 2.0,
                     labSize = 3.0,
                     labCol = 'black',
                     labFace = 'plain',
                     boxedLabels = TRUE,
                     colAlpha = 4/5,
                     legendPosition = 'right',
                     legendLabSize = 14,
                     legendIconSize = 4.0,
                     drawConnectors = TRUE,
                     widthConnectors = 1.0,
                     colConnectors = 'black')

