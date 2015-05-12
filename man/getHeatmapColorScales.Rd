\name{getHeatmapColorScales}
\alias{getHeatmapColorScales}
\title{
    Get caOmicsV Heatmap Color Scales
}
\description{
    Generate a color map for heatmap color scales
}
\usage{
    getHeatmapColorScales(colorType)
}
\arguments{
    \item{colorType}{
    character vector, one of "BlueWhiteRed", "GreenWhiteRed", 
    "GreenYellowRed", "GreenBlackRed", "YellowToRed", and "Black".
    }
}
\value{
    An RGB color matrix with dimension of 255 by 2.
}
\author{
    Henry Zhang
}
\examples{
    colorMap <- getHeatmapColorScales("BlueWhiteRed")
}
\keyword{ methods }
