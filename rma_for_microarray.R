#### Step1. install packages----Demo (GSE1297)####

# install the core bioconductor packages, if not already installed
source("http://bioconductor.org/biocLite.R")
biocLite()

# install additional bioconductor libraries, if not already installed
biocLite("GEOquery")  # Get data from NCBI Gene Expression Omnibus (GEO)
biocLite("affy")  # Methods for Affymetrix Oligonucleotide Arrays
biocLite("hgu133a.db", type = "source")  # GSE1297: Platform_title = [HG-U133A]
biocLite("hgu133acdf")

#### Step2. download the CEL file ####

library(GEOquery)

# Set working directory for download
setwd("D:/workspace/R/rawdata")

# Download the CEL file (by GSE - Geo series id), may take very long time
getGEOSuppFiles("GSE10961")
# this command does not work for this chip, so I need to download it
# directly Platform_title = [HG-U133A] Affymetrix Human Genome U133A Array
# Platform information can be found in SOFT file

# Unpack the CEL files
setwd("D:/workspace/R/rawdata/GSE40367")

untar("GSE40367_RAW.tar", exdir = "data")
cels = list.files("./data/", pattern = "CEL", ignore.case = TRUE)
# sometiles, it is 'CEL', you need to check it first
sapply(paste("data", cels, sep = "/"), gunzip)

cels = list.files("./data/", pattern = "CEL", ignore.case = TRUE)

#cels = list.files("./GPL93/", pattern = "CEL")
# sometiles, it is 'CEL', you need to check it first

#### Step3. Perform RMA normalization ####

library(affy)
library(hgu133a.db)
library(hgu133acdf)
###
library(u133x3p.db)
library(u133x3pcdf)

# Set working directory for normalization
setwd("D:/Work/Experssion/Breast cancer/GSE46141/data/")
###common method###
#raw.data = ReadAffy(verbose = FALSE, filenames = cels, cdfname = "hgu133acdf")
raw.data = ReadAffy(verbose = FALSE, filenames = cels)

###oligo method###
# load the oligo library
library(oligo)

# Read in the CEL files in the directory
celFiles <- list.celfiles()
raw.data <- read.celfiles(celFiles, phenoData = phenoData(gset[[1]]))
###read cel method end###

# perform RMA normalization (log2)
data.rma.norm = rma(raw.data)

# Get the expression estimates for each array
data.exprs = exprs(data.rma.norm)

# Take a look at the result (first 5 rows and 5 columes)
data.exprs[1:5, 1:5]

# Write RMA-normalized, mapped data to file
write.table(data.exprs, file = "GSE43837.txt", quote = FALSE, sep = "\t")


#### Step4. Annotation ####

tt = cbind(row.names(data.exprs), data.exprs)
###
library(stringr)
columnnames <- str_extract(colnames(data.exprs),pattern = "GSM\\d+")
###
colnames(tt) <- c("ProbID", columnnames)
colnames(tt) <- c("ProbID", sub(".CEL", "", colnames(data.exprs), ignore.case = TRUE))
#colnames(tt) = c("ProbID", sub("_HuGene-1_0-st-v1_.*.CEL", "", colnames(data.exprs), ignore.case = TRUE))
rownames(tt) = NULL
tt[1:5, 1:5]

require(RCurl)
myURL <- getURL("https://dl.dropboxusercontent.com/u/8272421/geo/HGU133A.na33.txt",
                ssl.verifypeer = FALSE)
annot <- read.table(textConnection(myURL), header = TRUE, sep = "\t")
head(annot)

####
library(hgu133plus2.db) # GPL570
columns(hgu133plus2.db)
ls("package:hgu133plus2.db")
annot <- toTable(hgu133plus2ENTREZID)
####
library(u133x3p.db) # GPL1352
annot <- toTable(u133x3pENTREZID)
####
library(hgu133b.db) # GPL97
annot <- toTable(hgu133bENTREZID)
####
annot <- read.csv("D:/workspace/R/rawdata/GSE90814/GPL5175.csv")
####

# probe sets were mapped to Entrez Gene IDs.
# comb=merge(annot,tt,by.x='ProbeSetID',by.y='ProbID',all.y=TRUE)
#comb = merge(annot, tt, by.x = "ID", by.y = "ProbID")
comb <- merge(annot, tt, by.x = "probe_id", by.y = "ProbID")
comb[1:5, 1:5]

write.table(comb, file = "comb2.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# If multiple probe sets corresponded to the same gene, then the expression
# values of these probe sets were averaged.
#comb2 <- subset(comb, select = -c(ProbeSetID))
comb2 <- subset(comb, select = -c(probe_id))
comb2 <- data.frame(lapply(comb2, as.character), stringsAsFactors = FALSE)
comb2 <- data.frame(lapply(comb2, as.numeric), stringsAsFactors = FALSE)
out <- aggregate(. ~ gene_id, data = comb2, mean)

# Format values to 5 decimal places
out <- format(out, digits = 5)
out[1:5, 1:5]

setwd("D:/workspace/R/rawdata/GSE40367/")

#write.table(out, file = "GSE43837.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(out, file = "GSE40367.csv", quote = FALSE, sep = ",", row.names = FALSE)

#### Step5. Session information ####

sessionInfo()

#### StepX1. extract platform information which contain probe-symbol-entrezid####
gset <- getGEO("GPL15207", destdir = './')
probetoentrezid <- Table(gset)[,c(1,15,19)]
annot <- probetoentrezid[!grepl('///|---', probetoentrezid$`Entrez Gene`), c("ID", "Entrez Gene")]

#### StepX2 ####
rm(list = ls())

#### Make cdf package by cdf file ####

biocLite("makecdfenv") #install makecdfenv package

cdfpath <- './data/'
pkgpath <- 'C:/Users/dxs/Documents/R/win-library/3.5'
pkgname <- 'hursta2a520709cdf'
make.cdf.package(filename = "hursta2a520709.gz", cdf.path = cdfpath, package.path = pkgpath,
                 packagename = pkgname, compress = TRUE, species = 'Homo_sapiens')
#use windows DOS or Linux shell
#Win+R change dir to cdf package then
#input `R CMD INSTALL hursta2a520709cdf`
#run