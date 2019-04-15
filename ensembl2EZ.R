# ENSEMBL ID to Gene ID(ENTREZID)
# or to Gene Symbol
library(org.Hs.eg.db)

columns(org.Hs.eg.db)

k <- keys(org.Hs.eg.db,keytype = "ENSEMBL")
k <- keys(org.Hs.eg.db,keytype = "SYMBOL")
head(k,5)

#ees <- select(org.Hs.eg.db,keys=k,columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
es <- select(org.Hs.eg.db,keys=k,columns = c("ENTREZID"), keytype="SYMBOL")
ee <- select(org.Hs.eg.db,keys=k,columns = c("ENTREZID"), keytype="ENSEMBL")

dim(ees)
head(ees,5)

FPKM <- as.matrix(read.table('GSE118435', sep = '\t', header = TRUE))
#ID <- FPKM[,1]
#ID_list <- ees[match(ID,ees[,"SYMBOL"]),]

comb <- merge(ees, FPKM, by.x = "SYMBOL", by.y = "Ensemble_ID")
comb2 <- subset(comb, select = -c(SYMBOL))

comb <- merge(ee, FPKM, by.x = "ENSEMBL", by.y = "ENSEMBL_ID")
comb2 <- subset(comb, select = -c(ENSEMBL))

comb2 <- data.frame(lapply(comb2, as.character), stringsAsFactors = FALSE)
comb2 <- data.frame(lapply(comb2, as.numeric), stringsAsFactors = FALSE)

out <- aggregate(. ~ ENTREZID, data = comb2, mean)

write.table(comb2, file = "GSE118435.txt", quote = FALSE, sep = "\t", row.names = FALSE)

colnames(comb2) <- c('ENTREZID', 'GSM1417182', 'GSM1417177', 'GSM1417184', 
                     'GSM1417181', 'GSM1417179', 'GSM1417183', 'GSM1417180', 'GSM1417178')
