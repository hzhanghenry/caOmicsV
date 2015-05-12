\name{plotBioMatrixSampleData}
\alias{plotBioMatrixSampleData}
\title{ 
    Plot Sample Data on caOmicsV bioMatrix Layout 
}
\description{
    Plot colored polygons on phenotype area of caOmicsV bioMatrix layout to 
    show sample information such as diagnosis, tissue type, ...  bioMatrix 
    layout and graphic device must be initialized first.
}
\usage{
    plotBioMatrixSampleData(rowNumber, areaName, fillColor=NA, 
        borderColor=NA, topAdjust=0, bottomAdjust=0)
}
\arguments{
    \item{rowNumber}{
    non-negative integer, number of the row where the data to be plotted
    }
    \item{areaName}{
    character vector, either "phenotype" or "omicsdata"
    }
    \item{fillColor}{
    character vector of color names or vector of R color specification
    }
    \item{borderColor}{
    character vector or a R colors specification for boarder color
    }
    \item{topAdjust}{
    non-negative numeric, height that will be reduced from top 
    }
    \item{bottomAdjust}{
    non-negative numeric, height that will be reduced from bottom 
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

    sampleColor <- c(rep("blue", 25), rep("red", 25))
    plotBioMatrixSampleData(rowNumber=1, sampleColor, areaName="phenotype") 
}
\keyword{ methods }
