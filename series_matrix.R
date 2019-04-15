##download GSExxxx series matrix
##extract sample information
library(GEOquery)
setwd("D:/Work/Sequencing/GSE50760_FPKM/")

dir.create("D:/workspace/R/rawdata/GSE50630")
setwd("D:/workspace/R/rawdata/GSE50630")

gset <- getGEO("GSE50760",destdir = './')
pdata <- phenoData(gset[[1]])
da <- as.data.frame(pdata@data)
##met_info <- as.matrix(da$`Tissue:ch1`)
met_info <- da[c('geo_accession','tumor site:ch1')]
met_info$`tumor site:ch1` <- gsub(", DMSO control","",met_info$`tumor site:ch1`)

#setwd("D:/workspace/R/rawdata/GSE50630")
write.table(met_info, file = "GSE107000_label.csv", quote = FALSE, sep = ",", row.names = FALSE)
##################################
exp <- exprs(gset[[1]])
fdata <- featureData(gset[[1]])
fd <- fdata@data
probe2Gene <- fd[c('ID', 'Entrez_Gene_ID')]
probetogene <- probe2Gene[grepl('', probe2Gene$Entrez_Gene_ID),]
write.table(probetogene, file = "GPL6947.csv", quote = FALSE, sep = ",", row.names = FALSE)
##################################
tt = cbind(row.names(exp), exp)
colnames(tt) <- c("ProbID", colnames(exp))
comb <- merge(probetogene, tt, by.x = "ID", by.y = "ProbID")
comb2 <- subset(comb, select = -c(ID))
comb2 <- data.frame(lapply(comb2, as.character), stringsAsFactors = FALSE)
comb2 <- data.frame(lapply(comb2, as.numeric), stringsAsFactors = FALSE)
out <- aggregate(. ~ Entrez_Gene_ID, data = comb2, mean)
out <- format(out, digits = 5)
out[1:5, 1:5]

write.table(out, file = "GSE50630.csv", quote = FALSE, sep = ",", row.names = FALSE)
