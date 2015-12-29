\name{setBioMatrixPlotParameters}
\alias{setBioMatrixPlotParameters}
\title{
    Set up Plot Parameters for caOmicsV bioMatrix Layout
}
\description{
    Put biomatrix plot parameters to CA_OMICS_ENV environment. All arguments 
    should be validated outside in advance. This function is for internal 
    use only.
}
\usage{
    setBioMatrixPlotParameters(numOfGenes, numOfSamples, numOfPhenotypes, 
            sampleHeight, sampleWidth, columnPadding, rowPadding, 
            geneNameWidth, remarkWidth, summaryWidth, sampleNameHeight, 
            legendHeight)
}
\arguments{
    \item{numOfGenes}{
    non-negative integer, total number of genes to be plotted
    }
    \item{numOfSamples}{
    non-negative integer, total number of samples to be plotted
    }
    \item{numOfPhenotypes}{
    non-negative integer, total number of phenotypes to be plotted
    }
    \item{sampleHeight}{
    non-negative numeric, height of polygon (sample) in inch on each row
    }
    \item{sampleWidth}{
    non-negative numeric, width of polygon (sample) in inch on each row
    }
    \item{columnPadding}{
    non-negative numeric, padding width in inch between two polygons (samples)
    }
    \item{rowPadding}{
    non-negative numeric, padding height in inch between two rows
    }
    \item{geneNameWidth}{
    non-negative numeric, width of left text plot area in inch for gene names
    }
    \item{remarkWidth}{
    non-negative numeric,  width of plot area in inch for right side text plot 
    (second set of gene names) and/or summary data plot
    } 
    \item{summaryWidth}{
    non-negative numeric,  width of plot area in inch on right side of data
    plot area for summary info plot
    } 
    \item{sampleNameHeight}{
    non-negative numeric, height of plot area in inch for sample names
    }
    \item{legendHeight}{ non-negative numeric, height of legend area }
}
\value{
    None
}
\author{
    Henry Zhang
}
\examples{
    #  Internal use only
}
\keyword{ methods }
