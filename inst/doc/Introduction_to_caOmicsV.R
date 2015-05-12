### R code from vignette source 'Introduction_to_caOmicsV.Rnw'

###################################################
### code chunk number 1: caOmicsVBioMatrixLayoutDemo
###################################################
library(caOmicsV)
data(caOmicsV.biomatrix.eset)
png("caOmicsV.bioMatrix.Layout.Plot.png", height=8, width=11,
    units="in", res=300, type="cairo")
plotBioMatrix(caOmicsV.biomatrix.eset, summaryType="text")
bioMatrixLegend(heatmapNames=c("RNASeq", "miRNASeq"),
    categoryNames=c("Methyl H", "Methyl L"),
    binaryNames=c("CN LOSS", "CN Gain"),
    heatmapMin=-3, heatmapMax=3, colorType="BlueWhiteRed")
dev.off()


###################################################
### code chunk number 2: caOmicsVbionetCircosLayoutDemo
###################################################
data(caOmicsV.bionet.eset)
png("caOmicsV.bioNetCircos.Layout.Plot.png", height=8, width=11,
    units="in", res=300, type="cairo")
plotBioNetCircos(caOmicsV.bionet.eset)
dataNames <- c("Tissue Type", "RNASeq", "miRNASeq", "Methylation", "CNV")
bioNetLegend(dataNames, heatmapMin=-3, heatmapMax=3)
dev.off()


###################################################
### code chunk number 3: caOmicsV_bionet_eset
###################################################
data(caOmicsV.biomatrix.eset)
names(caOmicsV.biomatrix.eset)

data(caOmicsV.bionet.eset)
names(caOmicsV.bionet.eset)


###################################################
### code chunk number 4: bioMatrixDemoData
###################################################
library(caOmicsV)
data(caOmicsV.biomatrix.eset)
eset <- caOmicsV.biomatrix.eset
names(eset)


###################################################
### code chunk number 5: InitializeBioMatrix
###################################################
numOfGenes <- length(eset$geneNames);
numOfSamples <- length(eset$sampleNames); 
numOfPhenotypes <- nrow(eset$sampleInfo)-1;

numOfHeatmap <- length(eset$heatmapData);
numOfSummary <- length(eset$summaryData);
phenotypes   <- rownames(eset$sampleInfo)[-1];

sampleHeight <- 0.4;
sampleWidth <- 0.1; 
samplePadding <- 0.025;
geneNameWidth <- 1;
sampleNameHeight <- 2.5;
remarkWidth <- numOfHeatmap; 
rowPadding <- 0.1;

initializeBioMatrixPlot(numOfGenes, numOfSamples, numOfPhenotypes, 
    sampleHeight, sampleWidth, samplePadding, rowPadding, 
    geneNameWidth, remarkWidth, sampleNameHeight)
caOmicsVColors <- getCaOmicsVColors()

png("caOmicsVbioMatrixLayoutDemo.png", height=8, width=12, 
        unit="in", res=300, type="cairo")
par(cex=0.75)
showBioMatrixPlotLayout(eset$geneNames,eset$sampleNames, phenotypes)


###################################################
### code chunk number 6: PlotBioMatrixPhenotype
###################################################
head(eset$sampleInfo)[,1:3]
rowIndex <- 2;

sampleGroup <- as.character(eset$sampleInfo[rowIndex,])
sampleTypes <- unique(sampleGroup)
sampleColors <- rep("blue", length(sampleGroup));
sampleColors[grep("Tumor", sampleGroup)] <- "red"

rowNumber <- 1
areaName <- "phenotype"
plotBioMatrixSampleData(rowNumber, areaName, sampleColors);

geneLabelX <- getBioMatrixGeneLabelWidth()
maxAreaX <- getBioMatrixDataAreaWidth()
legendH <- getBioMatrixLegendHeight()
plotAreaH <- getBioMatrixPlotAreaHeigth()
sampleH<- getBioMatrixSampleHeight()

sampleLegendX <- geneLabelX + maxAreaX
sampleLegendY <- plotAreaH + legendH - length(sampleTypes)*sampleH
colors <- c("blue", "red")
legend(sampleLegendX, sampleLegendY, legend=sampleTypes, 
    fill=colors,  bty="n", xjust=0)


###################################################
### code chunk number 7: PlotBioMatrixHeatmap
###################################################
heatmapData <- as.matrix(eset$heatmapData[[1]][,]);
plotBioMatrixHeatmap(heatmapData, maxValue=3, minValue=-3)

heatmapData <- as.matrix(eset$heatmapData[[2]][,])
plotBioMatrixHeatmap(heatmapData, topAdjust=sampleH/2,  
    maxValue=3, minValue=-3);

secondNames <- as.character(eset$secondGeneNames)
textColors <- rep(caOmicsVColors[3], length(secondNames));
plotBioMatrixRowNames(secondNames, "omicsData", textColors, 
    side="right", skipPlotColumns=0);


###################################################
### code chunk number 8: PlotBioMatrixCategoryData
###################################################
categoryData <- eset$categoryData[[1]]
totalCategory <- length(unique(as.numeric(eset$categoryData[[1]])))

plotColors <- rev(getCaOmicsVColors())
plotBioMatrixCategoryData(categoryData, areaName="omicsData", 
    sampleColors=plotColors[1:totalCategory])


###################################################
### code chunk number 9: PlotBioMatrixBinaryData
###################################################
binaryData <- eset$binaryData[[1]];
plotBioMatrixBinaryData(binaryData, sampleColor=caOmicsVColors[4]);

binaryData <- eset$binaryData[[2]];
plotBioMatrixBinaryData(binaryData, sampleColor=caOmicsVColors[3])


###################################################
### code chunk number 10: PlotBioMatrixSummaryData
###################################################
summaryData  <- eset$summaryInfo[[1]][, 2];
summaryTitle <- colnames(eset$summaryInfo[[1]])[2];

remarkWidth <- getBioMatrixRemarkWidth();
sampleWidth <- getBioMatrixSampleWidth();
col2skip <- remarkWidth/2/sampleWidth + 2;

plotBioMatrixRowNames(summaryTitle, areaName="phenotype", 
    colors="black", side="right", skipPlotColumns=col2skip);

plotBioMatrixRowNames(summaryData, "omicsData", 
    colors=caOmicsVColors[3], side="right", 
    skipPlotColumns=col2skip)


###################################################
### code chunk number 11: AddBioMatrixLegend
###################################################
bioMatrixLegend(heatmapNames=c("RNASeq", "miRNASeq"), 
    categoryNames=c("Methyl H", "Methyl L"), 
    binaryNames=c("CN LOSS", "CN Gain"),   
    heatmapMin=-3, heatmapMax=3, 
    colorType="BlueWhiteRed")
dev.off()


###################################################
### code chunk number 12: BioNetCircosDemoData
###################################################
library(caOmicsV)

data(caOmicsV.bionet.eset)
eset <- caOmicsV.bionet.eset

sampleNames  <- eset$sampleNames
geneNames    <- eset$geneNames
numOfSamples <- length(sampleNames)

numOfSampleInfo <- nrow(eset$sampleInfo) - 1
numOfSummary <- ifelse(eset$summaryByRow, 0, col(eset$summaryInfo)-1)
numOfHeatmap <- length(eset$heatmapData)
numOfCategory <- length(eset$categoryData)
numOfBinary <- length(eset$binaryData)

expr <- eset$heatmapData[[1]]
bioNet <- bc3net(expr) 


###################################################
### code chunk number 13: InitializeBioNetCircos
###################################################
widthOfSample    <- 100
widthBetweenNode <- 3
lengthOfRadius   <- 10

dataNum <- sum(numOfSampleInfo, numOfSummary, numOfHeatmap, 
    numOfCategory, numOfBinary)
trackheight <- 1.5
widthOfPlotArea  <- dataNum*2*trackheight

initializeBioNetCircos(bioNet, numOfSamples, widthOfSample, 
    lengthOfRadius, widthBetweenNode, widthOfPlotArea)
caOmicsVColors <- getCaOmicsVColors()
supportedType <- getCaOmicsVPlotTypes()

par(cex=0.75)
showBioNetNodesLayout()


###################################################
### code chunk number 14: BioNetCircosNodeLayout
###################################################
par(cex=0.6)
onTop <- c(14, 15, 16, 9, 7, 20, 8, 24, 10, 25)
labelBioNetNodeNames(nodeList=onTop,labelColor="blue", 
    labelLocation="top", labelOffset = 0.7)

onBottom <- c(26, 22, 23, 18, 19, 3, 5)
labelBioNetNodeNames(nodeList=onBottom,labelColor="black", 
    labelLocation="bottom", labelOffset = 0.7)

onLeft <- c(2, 11, 21, 17)
labelBioNetNodeNames(nodeList=onLeft,labelColor="red", 
    labelLocation="left", labelOffset = 0.7)

onRight <- c(13, 12, 4, 1, 6)
labelBioNetNodeNames(nodeList=onRight,labelColor="brown", 
    labelLocation="right", labelOffset = 0.7)


###################################################
### code chunk number 15: EraseBioNetCircosNodes
###################################################
eraseBioNetNode()


###################################################
### code chunk number 16: BioNetCircosBoundary
###################################################
inner <- lengthOfRadius/2
outer <- inner +  trackheight


###################################################
### code chunk number 17: PlotBioNetCircosPhenotype
###################################################
groupInfo <- as.character(eset$sampleInfo[2, ])
sampleColors <- rep("blue", numOfSamples);
sampleColors[grep("Tumor", groupInfo)] <- "red"

plotType=supportedType[1]
groupInfo <- matrix(groupInfo, nrow=1)
bioNetCircosPlot(dataValues=groupInfo, 
    plotType, outer, inner, sampleColors)

inner <- outer + 0.5
outer <- inner +  trackheight  


###################################################
### code chunk number 18: PlotBioNetCircosHeatmap
###################################################
exprData <- eset$heatmapData[[1]]
plotType <- supportedType[4] 
bioNetCircosPlot(exprData, plotType, outer, inner, 
    plotColors="BlueWhiteRed", maxValue=3, minValue=-3)

inner <- outer + 0.5 
outer <- inner +  trackheight


###################################################
### code chunk number 19: PlotBioNetCircosCategoryData
###################################################
categoryData <- eset$categoryData[[1]]
plotType <- supportedType[2];
bioNetCircosPlot(categoryData, plotType, 
    outer, inner, plotColors="red")
inner <- outer + 0.5
outer <- inner +  trackheight


###################################################
### code chunk number 20: PlotBioNetCircosBinaryData
###################################################
binaryData <- eset$binaryData[[1]]
plotType <- supportedType[3]
plotColors <- rep(caOmicsVColors[1], ncol(binaryData))
bioNetCircosPlot(binaryData, plotType, 
    outer, inner, plotColors)

inner <- outer + 0.5
outer <- inner +  trackheight 


###################################################
### code chunk number 21: LinkBioNetCircosSamples
###################################################
outer <- 2.5
bioNetGraph <- getBioNetGraph()
nodeIndex <- which(V(bioNetGraph)$name=="PLVAP")

fromSample <- 10 
toSample <- 50 
plotColors <- "red"
linkBioNetSamples(nodeIndex, fromSample, 
    toSample, outer, plotColors)

fromSample <- 40 
toSample <- 20 
plotColors <- "blue"
linkBioNetSamples(nodeIndex, fromSample, 
    toSample, outer, plotColors)


###################################################
### code chunk number 22: AddBioNetCircosLegend
###################################################
dataNames <- c("Tissue Type", "RNASeq", "Methylation", "CNV")
bioNetLegend(dataNames, heatmapMin=-3, heatmapMax=3)


###################################################
### code chunk number 23: sessionInfo
###################################################
sessionInfo()


