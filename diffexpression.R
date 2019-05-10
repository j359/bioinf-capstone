load("/Users/tesssenftle/Downloads/LUAD_mRNA.Rda")
load("/Users/tesssenftle/Downloads/LUAD_mRNA.subtypes.Rda")







library(plyr)
library(edgeR)



demrna <- t(finalmrna_filtered)
y <- DGEList(counts=demrna, group = cluster_outcome$group)

#design <- model.matrix(~group)
#fit <- glmQLFit(y, design)




#Filtering
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

#Normalization
y <- calcNormFactors(y)

#GLM
design <- model.matrix(~0+group, data = y$samples)
colnames(design) <- levels(y$samples$group)



y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmQLFit(y, design)



 qlf <- glmQLFTest(fit, contrast=c(-1,1,0))
 topTags(qlf)









#compares 1 and 2
lrt <- glmLRT(fit, contrast = c(-1,1,0,0))


c1and2 <- summary(de1 <- decideTestsDGE(lrt, p=0.05, adjust="BH"))

#compares 2 and 3
lrt <- glmLRT(fit, contrast = c(0,-1,1,0))
c2and3 <- summary(de1 <- decideTestsDGE(lrt, p=0.05, adjust="BH"))

#compares 3 and 4
lrt <- glmLRT(fit, contrast = c(0,0,-1,1))
c3and4 <- summary(de1 <- decideTestsDGE(lrt, p=0.05, adjust="BH"))

#compares 1 and 4
lrt <- glmLRT(fit, contrast = c(-1,0,0,1))
c1and4 <-summary(de1 <- decideTestsDGE(lrt, p=0.05, adjust="BH"))

#compares 2 and 4
lrt <- glmLRT(fit, contrast = c(0,-1,0,1))
c2and4 <- summary(de1 <- decideTestsDGE(lrt, p=0.05, adjust="BH"))

#compares 1 and 3
lrt <- glmLRT(fit, contrast = c(-1,0,1,0))

c1and3 <- summary(de1 <- decideTestsDGE(lrt, p=0.05, adjust="BH"))


detags1 <- rownames(y)[as.logical(de1)]

results1 = topTags(lrt, n = length(detags1))$table

nrow(results1) 

differentiated.genes1and3 <- rownames(results1)

#combine 

help(append)


combined1 <- append(differentiated.genes1and2, differentiated.genes1and3)

combined2 <- append(combined1, differentiated.genes1and4)

combined3 <- append(combined2, differentiated.genes2and4)

combined4 <- append(combined3, differentiated.genes2and3)

combined5 <- append(combined4, differentiated.genes3and4)



#write.csv(finalde, "/Users/tesssenftle/Downloads/finalde.csv")

finalde <- unique(combined5)

finalgenes <- read.csv(file = "/Users/tesssenftle/Downloads/finalde.csv")























