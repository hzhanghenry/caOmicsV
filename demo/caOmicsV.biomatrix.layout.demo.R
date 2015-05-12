#
#    Demo of caOmicsV bioMatrix Layout 
#
#    Last revised on March 6, 2015
#    ____________________________________________________________
#    <bioMatrixDemo><bioMatrixDemo><bioMatrixDemo><bioMatrixDemo>


library(caOmicsV)
data(caOmicsV.biomatrix.eset)

pdf("caOmicsV.bioMatrix.Layout.Plot.pdf", height=8, width=11)

plotBioMatrix(caOmicsV.biomatrix.eset, summaryType="text")
bioMatrixLegend(heatmapNames=c("RNASeq", "miRNASeq"), 
    categoryNames=c("Methyl H", "Methyl ZL"), 
    binaryNames=c("CN LOSS", "CN Gain"),   
    heatmapMin=-3, heatmapMax=3, colorType="BlueWhiteRed")

dev.off()