#
#    Demo of caOmicsV bioMatrix Layout 
#
#    Last revised on March 6, 2015
#    ____________________________________________________________
#    <bioMatrixDemo><bioMatrixDemo><bioMatrixDemo><bioMatrixDemo>


library(caOmicsV)
data(biomatrixPlotDemoData)

pdf("caOmicsV.bioMatrix.Layout.Plot.pdf", height=12, width=12)

plotBioMatrix(biomatrixPlotDemoData, summaryType="text")
bioMatrixLegend(heatmapNames=c("RNASeq", "miRNASeq"), 
    categoryNames=c("Methyl H", "Methyl L"), 
    binaryNames=c("CN LOSS", "CN Gain"),   
    heatmapMin=-3, heatmapMax=3, colorType="BlueWhiteRed")

dev.off()
