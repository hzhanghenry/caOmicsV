\name{plotHeatmapColorScale}
\alias{plotHeatmapColorScale}
\title{
    Plot Heatmap Color Scale for both caOmicsV bioMatrix and bioNewCircos 
    Layout
}
\description{
    Draw heatmap color scale for legends of bioMatrix plot and bioNetCircos 
    plot. Graphic device must be initialized.
}
\usage{
    plotHeatmapColorScale(coorX, coorY, colorType="BlueWhiteRed", 
        scaleWidth, scaleHeight, minValue, maxValue, direction="h")
}
\arguments{
    \item{coorX}{
    numeric, x coordinates for the top left location of color scale
    }
    \item{coorY}{
    numeric, y coordinates for the top left location of color scale
    }
    \item{scaleWidth}{ non-negative numeric, width of color scale }
    \item{scaleHeight}{ non-negative numeric, height of color scale }
    \item{colorType}{
    character vector, one of "BlueWhiteRed", "GreenWhiteRed", 
    "GreenYellowRed", "GreenBlackRed", "YellowToRed", and "BlackOnly"
    }
    \item{minValue}{
    numeric, the smallest value associated with the lowest coluor
    }
    \item{maxValue}{
    numeric, the biggest value associated with the highest coluor
    }
    \item{direction}{
    One character, either "h" for horizontal or "v" for vertical
    }
}
\value{
    None
}
\author{
    Henry Zhang
}
\examples{
    plot.new()
    plot.window(xlim=c(1,5), ylim=c(1,5))
    axis(1)
    axis(2)
    plotHeatmapColorScale(coorX=2, coorY=4, colorType="BlueWhiteRed", 
        scaleWidth=2, scaleHeight=0.5, minValue=-3, maxValue=3)
}
\keyword{ methods }
