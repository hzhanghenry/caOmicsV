#
#   Demo of caOmicsV bioNetCircos Layout 
#
#   Last revised on March 6, 2015
#    ________________________________________________________________________
#    <bioNetCircosDemo><bioNetCircosDemo><bioNetCircosDemo><bioNetCircosDemo>


library(caOmicsV)
data(bionetPlotDemoData)

pdf("caOmicsV.bioNetCircos.Layout.Plot.pdf", height=8, width=11)

plotBioNetCircos(bionetPlotDemoData)
dataNames <- c("Tissue Type", "RNASeq", "miRNASeq", "Methylation", "CNV")
bioNetLegend(dataNames, heatmapMin=-3, heatmapMax=3)

dev.off()
