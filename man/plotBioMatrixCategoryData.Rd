\name{plotBioMatrixCategoryData}
\alias{plotBioMatrixCategoryData}
\title{
    Plot Category Data on caOmicsV bioMatrix Layout
}
\description{
    Draw rectangle outline for one or more row(s) of samples to represent 
    categorical values. This function highlights all samples on each row. 
    bioMatrix layout and graphic device must be initialized first.
}
\usage{
    plotBioMatrixCategoryData(categoryData, 
            areaName=c("omicsData", "phenotype"), 
            sampleColors=palette(), 
            lineWidth=1, skipPlotRow=0)
}
\arguments{
    \item{categoryData}{
    vector or matrix of categorical values, such as 'High', 'low', and 'No'
    }
    \item{areaName}{
    character vector, either "omicsdata" or "phenotype"
    }
    \item{sampleColors}{
    character vector for color names or vector of R color specification
    }
    \item{lineWidth}{
    non-negative integer, graphic parameter for lwd (line width), default 1
    }
    \item{skipPlotRow}{
    non-negative integer, total rows on plot area that should be skipped when 
    only one row to plot, default 0. 
    }
}
\value{
    None
}
\author{
    Henry Zhang
}
\examples{
    initializeBioMatrixPlot(numOfGenes=1, numOfSamples=50)
    showBioMatrixPlotLayout("Gene", paste("Sample", 1:50),  "Diagnosis")

    categoryData <- matrix(c(rep(1, 15), rep(0, 20), rep(1, 15)), nrow=1)
    plotBioMatrixCategoryData(categoryData, areaName="omicsData") 
}
\keyword{ methods }
