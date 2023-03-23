library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-BLCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query,
            method = "api",
            files.per.chunk = 100) 
expression_data <- GDCprepare(query = query)
library(DT)
library(dplyr)
library(SummarizedExperiment)
count_data <- assay(expression_data)  ####下载完成，等待DeSeq2标准化
class(count_data)
count_data<- as.data.frame(count_data)
save(count_data,file="blca_mRNA_raw.Rdata")
count_data<- count_data+1
count_data <- apply(count_data, c(1,2), log2)
expr <- as.data.frame(count_data)
expr<- as.data.frame(t(expr))
expr<- cbind(rownames(expr),expr)
colnames(expr)[1] <- "ID"
expr$group <- substr(expr$ID,14,16)
expr <- subset(expr,expr$group=="01A")
expr$ID<- substr(expr$ID,1,12)
expr<- expr[!duplicated(expr$ID),]
###############
expr<- expr[,-56406]
expr<- as.data.frame(t(expr))
colnames(expr) <-expr[1,]
#############
exp<- data.table::fread("tcgaRisk.txt")
exp<- exp[,-(2:13)]
colnames(exp) <- c("ID","group")
expr<- as.data.frame(t(expr))
ee <- inner_join(exp,expr,by="ID")
ee<- as.data.frame(t(ee))
group<- as.data.frame(t(ee[1:2,]))
###
colnames(ee) <- ee[1,]
ee<- ee[-(1:2),]
###
group <- group$group 
group <- factor(group,levels = c("high","low"),ordered = F)
design <- model.matrix(~group)
colnames(design) <- levels(group)
design
library(limma)
fit <- lmFit(ee,design)
fit2 <- eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
diffLab <- subset(allDiff,abs(logFC) >1 & adj.P.Val < 0.05)
###
