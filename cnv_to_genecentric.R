#Preprocessing CNV data
library(reshape)
library(tidyr)
library(DNAcopy)
library(CNTools)

GDC_aliquot <- cnv.expression$GDC_Aliquot
Segment_Mean <- cnv.expression$Segment_Mean
Samples <- cnv.expression$Sample
Markers <- 1:191096
df <- cbind(Chromosome, Start, End, Segment_Mean, Samples)
df <- as.data.frame(df, stringsAsFactors=FALSE)
cnv_df <- spread(df, Samples, Segment_Mean)
cnv_df <- as.data.frame(lapply(cnv_df,as.numeric)) 
Start <- cnv_df$Start
Chromosome <- cnv_df$Chromosome
cnv_df <- cnv_df[,-c(1)]
sample_ids <- colnames(cnv_df)
cna.object <- CNA(cnv_df, Chromosome, Start, data.type="logratio",
                  sampleid = sample_ids, presorted = FALSE)
segout <- segment(cna.object)
sampleData <- segout[["output"]]
cnseg <- CNSeg(sampleData[which(is.element(sampleData[, "ID"], sample(unique(sampleData[, "ID"])))), ])
rdseg <- getRS(cnseg, by = "region", imput = FALSE, XY = FALSE, what = "mean")
data("geneInfo")
geneInfo <- geneInfo[sample(1:nrow(geneInfo)), ]
rdByGene <- getRS(cnseg, by = "gene", imput = FALSE, XY = FALSE, geneMap = geneInfo, what = "median")
reducedseg <- rs(rdByGene)
save(reducedseg, file = "cnv_data.Rda")
save(reducedseg, file = "cnv_data.csv")
