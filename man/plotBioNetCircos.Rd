\name{plotBioNetCircos}
\alias{plotBioNetCircos}
\title{
    Default Plot Method for caOmicsV bioNetCircos Layout
}
\description{
    A sample way to plot sample information and genomic data with default 
    settings. This method is mainly for demo purpose. In most cases, users 
    need make changes based on this display, specifically the igraph layout.
}
\usage{
    plotBioNetCircos(dataSet, graph=NULL, heatmapMax=NULL, heatmapMin=NULL, 
            heatmapColor="BlueWhiteRed")
}
\arguments{
    \item{dataSet}{
    an object returned from method of getESet() that contains all data for 
    plot
    }
    \item{graph}{
    an igraph object generated with one plot dataset, e.g, gene expression 
    value
    }
    \item{heatmapColor}{
    character vector, one of "BlueWhiteRed", "GreenWhiteRed", 
    "GreenYellowRed", "GreenBlackRed", "YellowToRed", and "Black".
    }
    \item{heatmapMax}{
    numeric, maximum value for heatmap plot, set to NULL to use the biggest 
    value in input data set.
    }
    \item{heatmapMin}{
    numeric, minimum value for heatmap plot, set to NULL to use the smallest 
    value in input data set.
    }
}
\value{
    None
}
\author{
    Henry Zhang
}
\examples{
    data(bionetPlotDemoData)
    plotBioNetCircos(bionetPlotDemoData)
}
\keyword{ methods }
