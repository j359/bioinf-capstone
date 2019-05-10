library(TCGAbiolinks)
library(SummarizedExperiment)

#clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
###################################################################
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Clinical",
                  file.type = "xml")
GDCdownload(query)
prepared.data <- GDCprepare_clinic(query, clinical.info = "patient")
#matrix = assay(prepared.data)

LUAD_CLINICAL_DATA<-prepared.data
save(LUAD_CLINICAL_DATA,file="LUAD_CLINICAL_DATA.Rda")

samples <- LUAD_CLINICAL_DATA$bcr_patient_barcode

#################################################################################################
mRNA.query <- GDCquery(project =  "TCGA-LUAD",
                       data.category = "Gene expression",
                       data.type = "Gene expression quantification",
                       experimental.strategy = "RNA-Seq",
                       file.type = "results",
                       barcode = samples,
                       legacy = TRUE)
GDCdownload(mRNA.query, method = 'api') #GDCdownload(mRNA.query)
mRNA.expression <- GDCprepare(mRNA.query)
mRNA = assay(mRNA.expression)

save(mRNA,file="LUAD_mRNA.Rda")

clinic = colData(mRNA.expression)
colnames(clinic)
subtype.of.samples = clinic[,c("sample", "patient", "barcode", "subtype_expression_subtype")]
save(subtype.of.samples, file = "LUAD_mRNA.subtypes.Rda")
#######################################################################################
miRNA.query <- GDCquery(project =  "TCGA-LUAD",
                        #experimental.strategy = "miRNA-Seq",
                        data.category = "Gene expression", 
                        barcode = samples,
                        data.type = "miRNA gene quantification",
                        platform = "Illumina HiSeq",
                        file.type = "hg19.mirbase20",
                        legacy = TRUE)
GDCdownload(miRNA.query, method = 'api') #GDCdownload(mRNA.query)
miRNA.expression <- GDCprepare(miRNA.query)
miRNA = assay(miRNA.expression)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("illuminaHumanv4.db")
library(illuminaHumanv4.db)

miRNAprobes <- miRNA.expression$miRNA_ID
miRNAnames <- data.frame(Gene=unlist(mget(x = miRNAprobes,envir = illuminaHumanv4SYMBOL)))

save(mRNA,file="LUAD_mRNA.Rda")

clinic = colData(mRNA.expression)
colnames(clinic)

########################################################################################
methyl.query <- GDCquery(project =  "TCGA-LUAD",
                       data.category = "DNA methylation",
                       #data.type = "Gene expression quantification",
                       #experimental.strategy = "RNA-Seq",
                       #file.type = "results",
                       platform = "Illumina Human Methylation 450",
                       barcode = samples,
                       legacy = TRUE)
GDCdownload(methyl.query, method = 'api') #GDCdownload(mRNA.query)
methyl.expression <- GDCprepare(methyl.query)
methyl = assay(methyl.expression)

save(methyl,file="LUAD_methyl.Rda")
write.csv(methyl,file="LUAD_methyl.csv")

clinic = colData(mRNA.expression)
colnames(clinic)

######################################################################################
cnv.query <- GDCquery(project =  "TCGA-LUAD",
                         data.category = "Copy Number Variation",
                         data.type = "Masked Copy Number Segment",
                         #experimental.strategy = "Genotyping Array",
                         #file.type = "hg19.seg",
                         platform = "Affymetrix SNP 6.0",
                         workflow.type = "DNAcopy",
                         barcode = samples)
GDCdownload(cnv.query, method = 'api') #GDCdownload(mRNA.query)
cnv.expression <- GDCprepare(cnv.query)
cnv = assay(cnv.expression)

save(cnv,file="LUAD_cnv.Rda")

