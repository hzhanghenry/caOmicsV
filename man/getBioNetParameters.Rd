\name{getBioNetParameters}
\alias{getBioNetBasePositions
    getBioNetGraph
    getBioNetNodePaddingScale
    getBioNetNodeParameters
    getBioNetNodePlotAreaBoundary
    getBioNetNodeRadius
    getBioNetPlotAreaWidth
    getBioNetPlotSampleWidth
    getBioNetPlotTotalSample
}
\title{
    Methods to Get caOmicsV BioNetCircos Plot Parameters
}
\description{
    Methods for retrieving plot parameters of caOmicsV bioNetCircos layout 
    stored in caOmicsV environment. The bioNetCircos layout must be 
    initialized first.
}
\usage{
    getBioNetBasePositions()
    getBioNetGraph()
    getBioNetNodePaddingScale()
    getBioNetNodeParameters()
    getBioNetNodePlotAreaBoundary()
    getBioNetNodeRadius()
    getBioNetPlotAreaWidth()
    getBioNetPlotSampleWidth()
    getBioNetPlotTotalSample()
}
\value{
    getBioNetBasePositions() returns numeric matrix containing x and y 
    coordinates of points on a circular line with radius of 1 and degrees of 
    text rotation on each point.  

    getBioNetGraph() returns an igraph object which representing a biological 
    network built with user's inputs.

    getBioNetNodePaddingScale() returns a numeric value for padding between 
    two nodes.  

    getBioNetNodeParameters() returns a list containing totalSamples, 
    sampleWidth, nodeRadius, nodePadding, plotAreaWidth, inner, and outer 
    boundary of plot area.

    getBioNetNodePlotAreaBoundary() returns a numeric vector for outer and 
    inner boundary of a node on caOmicsV bioNetCircos layout.  

    getBioNetNodeRadius() returns the numeric value for radius of a node on 
    caOmicsV bioNetCircos layout. 

    getBioNetPlotAreaWidth() returns the width of bioNetCircos layout plot 
    area.  

    getBioNetPlotSampleWidth() returns total points a sample will need on a 
    circular track.  

    getBioNetPlotTotalSample() returns the total number of samples to be 
    plotted on each node.  
}
\author{
    Henry Zhang
}
\examples{
    data(bionetPlotDemoData)
    expr <- bionetPlotDemoData$heatmapData[[1]]
    bioNet <- bc3net(expr) 
    initializeBioNetCircos(bioNet)

    positions   <- getBioNetBasePositions()
    bioGraph    <- getBioNetGraph()
    paddScale   <- getBioNetNodePaddingScale()
    nodeParams  <- getBioNetNodeParameters()
    plotBound   <- getBioNetNodePlotAreaBoundary()
    nodeRadius  <- getBioNetNodeRadius()
    areaWidth   <- getBioNetPlotAreaWidth()
    sampleWidth <- getBioNetPlotSampleWidth()
    numOfSample <- getBioNetPlotTotalSample()
}
\keyword{ methods }
