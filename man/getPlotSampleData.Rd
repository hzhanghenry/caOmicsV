\name{getPlotSampleData}
\alias{getPlotSampleData}
\title{
    Extract Subset of Sample Information
}
\description{
    Extract required rows and columns from a sample dataset
}
\usage{
    getPlotSampleData(sampleData, sampleNames)
}
\arguments{
    \item{sampleData}{
    Data frame with rows for samples and columns for features. Column 1 must 
    be sample names.
    }
    \item{sampleNames}{
    character vector, names of samples to select.
    }
}
\value{
    A data frame with subset of sample data and with the row order same as 
    sampleNames.
    }
\author{
    Henry Zhang
}
\examples{
    data(sampleDemoData)
    sampleNames <- as.character(sampleDemoData[10:40,1])
    sampleInfo <- getPlotSampleData(sampleDemoData, sampleNames) 
}
\keyword{ methods }
