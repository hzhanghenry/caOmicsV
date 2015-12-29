\name{showBioMatrixPlotLayout}
\alias{showBioMatrixPlotLayout}
\title{
    Display caOmicsV bioMatrix Layout
}
\description{
    Display caOmicsV bioMatrix layout with default information for each plot 
    area for purpose of optimizing the layout. A bioMatrix layout must be 
    initialized first.
}
\usage{
    showBioMatrixPlotLayout(geneNames, sampleNames, phenotypes, 
            secondGeneNames,  sampleColors=NULL, geneColors=NULL, 
            phenoColors=NULL)
}
\arguments{
    \item{geneNames}{
    character vector, gene names shown at left side of omics data plot area
    }
    \item{sampleNames}{
    character vector, sample names shown on the top of phenotype plot area
    }
    \item{phenotypes}{
    character vector, phenotype names listed on left side of phenotype plot 
    area
    }
    \item{secondGeneNames}{
    character vector, second gene names shown at right side of omics data area
    }
    \item{sampleColors}{
    character vector for color name(s) of R color specification for samples
    }
    \item{geneColors}{
    character vector of color name or R color specification for genes
    }
    \item{phenoColors}{
    character vector of color name(s) or R color specification for phenotypes
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
    showBioMatrixPlotLayout(geneNames="Gene", 
            sampleNames=paste("Sample", 1:50),  
            phenotypes="Diagnosis", secondGeneNames="miRNA")
}
\keyword{ methods }
