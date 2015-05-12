\name{plotBioMatrixRowNames}
\alias{plotBioMatrixRowNames}
\title{
    Plot Row Names on caOmicsV bioMatrix Layout
}
\description{
    Plot row names on the left or right side of biomatrix plot area. Phenotype 
    names, gene names, and remark notes are all plotted with this function. 
    bioMatrix layout and graphic device must be initialized first.
}
\usage{
    plotBioMatrixRowNames(geneNames, areaName, colors, side="left", 
            skipPlotRows=0, skipPlotColumns=0)
}
\arguments{
    \item{geneNames}{
    character vector, row names to be plotted
    }
    \item{areaName}{
    character vector, either "omicsData" or "phenotype"
    }
    \item{colors}{
    character vector of color names or R color specification
    }
    \item{side}{
    character vector, either "left" or "right"
    }
    \item{skipPlotRows}{
    non-negative integer, total rows on plot area that should be skipped from 
    the first row, default 0
    }
    \item{skipPlotColumns}{
    non-negative integer, columns (sampleWidth) will be skipped when plotting 
    items on remark area
    }
}
\value{
    None
}
\author{
    Henry Zhang
}
\examples{
    initializeBioMatrixPlot(numOfGenes=20, numOfSamples=50)

    geneNames  <- paste0("gene_", 1:20)
    miRNANames <- paste0("miRNA_", 1:20)

    showBioMatrixPlotLayout(geneNames, paste("Sample", 1:50),  "Diagnosis")
    plotBioMatrixRowNames(miRNANames, areaName="omicsData", colors="blue", 
                side="right")
}
\keyword{ methods }
