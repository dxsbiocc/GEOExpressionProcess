####使用R库oligo分析获取基因水平以及外显子水平的信号

library(oligo)
celFiles<-list.celfiles()#读取当前目录下所有的cel文件名
library(pd.huex.1.0.st.v2)#调取注释文件
affyRaw<-read.celfiles(celFiles)#读取CEL文件
dir.create("png")#创建一个目录，用于保存质量分析图片
for(i in 1:9){png(paste("png",i,"png",sep="."));image(affyRaw,i,col=gray((64:0)/64));dev.off();}#生成信号强度灰度图
pdf("png/hist.pdf")#生成柱状图
hist(affyRaw)
dev.off()
for(i in 1:9){png(paste("png/maplot",i,"png",sep=".")); MAplot(affyRaw,which=i,ylim=c(-1,1)); dev.off();}#生成原始信号MA图
png(paste("png/boxplot.png"))#生成原始信号箱线图
boxplot(affyRaw)
dev.off()
probesetSummaries<-rma(affyRaw, target="probeset")#计算外显子水平信号
geneSummaries<-rma(affyRaw,target="core")#计算基因水平信号
png("png/boxplot.rma.png")#生成RMA之后的箱线图
boxplot(geneSummaries,transfo=identity,main="After RMA")
dev.off()
pdf("png/hist.rma.pdf")#生成RMA之后的柱状图
hist(geneSummaries,transfo=identity,main="After RMA")
dev.off()
d.gene<-exprs(geneSummaries)#获得基因表达的数据
d.exon<-exprs(probesetSummaries)#获得外显子表达的数据

####使用R库xps分析获取基因水平以及外显子水平的信号

library(xps)
dir<-"/path/to/cel/data"#保存工作目录，因为xps在读取数据时可能会改变工作目录，所以需要随时将其改回来。
setwd(dir)
celfiles<-dir(pattern=".cel",ignore.case=T)#查找工作目录下所有的CEL文件
libdir <- "/path/to/library/HuEx-1_0-st-v2"#设置芯片注释库文件的目录
scheme.huex10stv2 <- import.exon.scheme("Scheme_HuEx10stv2",
                                        layoutfile=paste(libdir,"HuEx-1_0-st-v2.r2.clf",sep="/"),
                                        schemefile=paste(libdir,"HuEx-1_0-st-v2.r2.pgf",sep="/"),
                                        probeset=paste(libdir,"HuEx-1_0-st-v2.na31.hg19.probeset.csv",sep="/"),
                                        transcript=paste(libdir,"HuEx-1_0-st-v2.na31.hg19.transcript.csv",sep="/"))#读入注释库。xps采用了root，会将数据写入一个以.root结尾的文件当中，以便将来再次调用。这里就会在工作目录下生成一个名为Scheme_HuEx10stv2.root的文件。
celnames<-gsub("\\(HuEx-1_0-st-v2\\).CEL","",celfiles)#将colname做一个修改
data <- import.data(scheme.huex10stv2,filename="tmpdt_Data",filedir=dir,celdir=dir,celfiles=celfiles,celnames=celnames)#使用之前准备好的注释文件，读入CEL数据。
data.rma.ps<-rma(data,"HuExonRMAcore_ps",background="antigenomic",normalize=T,option="probeset",exonlevel="core")#在exon水平上做RMA，需要注意的是，xps不象其它的包那样对rma的结果做log2的转变，所以如果需要使用这里的值，将来还需要做一步log2 transform.
xps.rma.ps<-validData(data.rma.ps)#在上一步中，数据同时被写入一个名为HuExonRMAcore_ps.root的文件，方便将来调用。
data.rma.tc <- rma(data,"HuExonRMAcore_tc",background="antigenomic",normalize=T,option="transcript",exonlevel="core")#在gene水平上做RMA
xps.rma.tc<-validData(data.rma.tc)#数据被validData调用后才能正常操作，否则会显示为空数据
d.gene<-log2(xps.rma.tc)#获得基因表达的数据
d.exon<-log2(xps.rma.ps)#获得外显子表达的数据

####使用R库aroma.affymetrix分析获取基因水平以及外显子水平的信号
library(aroma.affymetrix)#调入库
verbose<-Arguments$getVerbose(-8,timestamp=TRUE)#设置返回信息时间格式
chipType<-"HuEx-1_0-st-v2"#设置芯片类型
cdf<-AffymetrixCdfFile$byChipType(chipType,tags="core,A20071112,EP")#设置芯片CDF文件
cs<-AffymetrixCelSet$byName("projectname_qiuworld",cdf=cdf)#设置CEL数据位置
bc <- RmaBackgroundCorrection(cs,tag="core")#设置RMA背噪参数
csBC <- process(bc,verbose=verbose)#DO RMA
qn <- QuantileNormalization(csBC, typesToUpdate="pm")#四分位标准化参数，其实是标准RMA的一个步骤
csN <- process(qn,verbose=verbose)#DO normalization
plmTr <- ExonRmaPlm(csN, mergeGroups=TRUE)#设置使用PLM模型获取基因水平表达数据
fit(plmTr,verbose=verbose)#获取基因水平表达数据
plmEx <- ExonRmaPlm(csN, mergeGroups=FALSE)#设置使用PLM模型获取外显子水平表达数据
fit(plmEx,verbose=verbose)#获取外显子水平表达数据
rs <-calculateResiduals(plmTr,verbose=verbose)#计算残基，用于绘制质量控制的NUSE图和RLE图
ae <- ArrayExplorer(rs)#读取数据并绘制图片
setColorMaps(ae,c("log2,log2neg,rainbow","log2,log2pos,rainbow"))
process(ae,interleaved="auto",verbose=verbose)
display(ae)
qamTr <- QualityAssessmentModel(plmTr)
plotNuse(qamTr)#NUSE图
plotRle(qamTr)#RLE图
cesTr <- getChipEffectSet(plmTr)#准备读入基因水平表达值
trFit <- extractDataFrame(cesTr, units=NULL, addNames=TRUE)#读入数据
cesEx <- getChipEffectSet(plmEx)#准备读入外显子表达值
exFit <- extractDataFrame(cesEx, units=NULL, addNames=TRUE)#读入数据
d.gene <- log2(trFit[,6:ncol(trFit)])#对数据做log2转换
row.names(d.gene)<-as.character(trFit[,1])#设置probeset_id为行号
d.exon <- log2(exFit[,6:ncol(trFit)])#对数据做log2转换
row.names(d.exon)<-as.character(exFit[,2])#设置probeset_id为行号

