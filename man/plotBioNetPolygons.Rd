\name{plotBioNetPolygons}
\alias{plotBioNetPolygons}
\title{
    Polygon Plot on caOmicsV bioNetCircos Layout
}
\description{
    Plot category data as polygons on node(s) of caOmicsV bioNetCircos layout. 
    bioNetCircos layout and graphic device  must be initialized first. Polygon 
    colours are converted from category value with build in colour series. Use 
    getcaOmicsVColors() and setcaOmicsVColors(colorList) to see or reset 
    customized colours.
}
\usage{
    plotBioNetPolygons(dataValues, outer, inner)
}
\arguments{
    \item{dataValues}{
    matrix of character or numeric for category data
    }
    \item{outer}{
    non-negative numeric, the outer boundary of plot area from node center
    }
    \item{inner}{
    non-negative numeric, the innner boundary of plot area from node center
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
    expr <- bionetPlotDemoData$heatmapData[[1]]
    bioNet <- bc3net(expr)

    initializeBioNetCircos(bioNet, totalSamples=60)
    showBioNetNodesLayout()
    dataValues <- matrix(rep(c(1:3), each=20), nrow=1)
    plotBioNetPolygons(dataValues, outer=3, inner=2)
}
\keyword{ methods }
