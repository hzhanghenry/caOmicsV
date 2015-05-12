\name{getPlotSummaryData}
\alias{getPlotSummaryData}
\title{
    Extract Summary Subset for Plotting
}
\description{
    Extract required rows and columns from a summary data set.
}
\usage{
    getPlotSummaryData(summaryData, sampleNames=NULL, geneNames=NULL)
}
\arguments{
    \item{summaryData}{
    a data frame with summary data for each gene (rows are for genes and 
    columns are summary values) or for each sample (rows are summary values 
    and columns are sample names)
    }
    \item{sampleNames}{
    character vector, names of samples/columns to be extracted. Set to NULL 
    when summary data is for genes.
    }
    \item{geneNames}{
    character vector, names of genes to be extracted. Set to NULL when summary 
    data is for samples
    }
}
\value{
    A data frame with subset of input data and with the same orders as 
    geneNames and sampleNames.
}
\author{
    Henry Zhang
}
\examples{
    data(RNASeqDemoData)
    exprValue <- as.matrix(RNASeqDemoData[, 2:ncol(RNASeqDemoData)])

    summaryData <- data.frame(as.character(RNASeqDemoData[,1]), 
        rowMeans(exprValue), log2(rowMeans(exprValue)))
    geneSymbols <- c("ECM1", "SLC26A6", "ADAMTS13", "FCN3", "CFP")
    sumByGene <- getPlotSummaryData(summaryData, geneNames=geneSymbols)

    summaryData <- rbind(colMeans(exprValue),  log2(colMeans(exprValue)))
    summaryData <- data.frame(c("sampleMean", "log2mean"), summaryData)
    colnames(summaryData) <- colnames(exprValue)
    sampleID <- colnames(exprValue)[c(2:21, 32:41)]
    sumBySample <- getPlotSummaryData(summaryData, sampleID)
}
\keyword{ methods }
