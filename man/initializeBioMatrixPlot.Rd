\name{initializeBioMatrixPlot}
\alias{initializeBioMatrixPlot}
\title{
    Set Up Parameters for caOmicsV bioMatrix Plot Layout
}
\description{
    Initialize parameters for caOmicsV bioMatrix layout to set up plot area 
    size and item sizes.  
} 
\usage{
    initializeBioMatrixPlot(numOfGenes=100, numOfSamples=100,   
        numOfPhenotypes=1, sampleHeight=0.4, sampleWidth=0.1,   
        columnPadding=0.025, rowPadding=0.1, geneNameWidth=1,   
        remarkWidth=1, summaryWidth=1, sampleNameHeight=1, 
        legendHeight=1)
}
\arguments{
    \item{numOfGenes}{
    non-negative numeric, total number of genes to be plotted
    }
    \item{numOfSamples}{
    non-negative numeric, total number of samples to be plotted
    }
    \item{numOfPhenotypes}{
    non-negative numeric, total number of phenotypes to be plotted
    }
    \item{sampleHeight}{
    non-negative numeric, height of rectangle area in inch for a sample plot
    }
    \item{sampleWidth}{
    non-negative numeric, width of rectangle area in inch for a sample plot
    }
    \item{columnPadding}{
    non-negative numeric, width of padding, in inch, between two rectangles 
    (samples)
    }
    \item{rowPadding}{
    non-negative numeric, height of padding in inch between two rows (genes)
    }
    \item{geneNameWidth}{
    non-negative numeric, width of plot area in inch for gene name plot
    }
    \item{remarkWidth}{
    non-negative numeric, width of plot area in inch for second set of gene 
    names and summary data.
    }
    \item{summaryWidth}{
    non-negative numeric, width of data plot area in inch to plot summary 
    info in text or bars.
    }
    \item{sampleNameHeight}{
    non-negative numeric, height of plot area in inch for sample names (sample 
    name are plotted vertically).
    }
    \item{legendHeight}{
    non-negative numeric, height of plot area in inch for legend.
    }
}
\value{
    None
}
\author{
    Henry Zhang
}
\examples{
    initializeBioMatrixPlot(numOfGenes=100, numOfSamples=100,
        numOfPhenotypes=1, sampleHeight=0.4, sampleWidth=0.1,
        columnPadding=0.025, rowPadding=0.1, geneNameWidth=1,
        remarkWidth=1, summaryWidth=1, sampleNameHeight=1, 
        legendHeight=1)
}
\keyword{ methods }
