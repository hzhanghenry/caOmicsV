\name{getPlotOmicsData}
\alias{getPlotOmicsData}
\title{
    Extract Subset from A Data Frame 
}
\description{
    Extract a subset from the input data for a set of samples and genes
}
\usage{
    getPlotOmicsData(omicsData, sampleNames, geneNames)
}
\arguments{
    \item{omicsData}{
    a data frame with all samples and all genes
    }
    \item{sampleNames}{
    character vector, names of samples to be extracted from dataset.
    }
    \item{geneNames}{
    character vector, names of genes to be extracted from dataset.
    }
}
\value{
    A data frame with subset of input data
}
\author{
    Henry Zhang
}
\examples{
    data(RNASeq)
    data(RNASeqDemoData)
    geneNames  <- as.character(RNASeqDemoData[,1])
    sampleNames <- colnames(RNASeqDemoData)[-1]
    plotData <- getPlotOmicsData(RNASeq, sampleNames, geneNames) 
}
\keyword{ methods }
