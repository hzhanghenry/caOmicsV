\name{bioMatrixLegend}
\alias{bioMatrixLegend}
\title{
    Plot Legend on caOmicsV bioMatrix Layout
}
\description{
    Draw legend including of heatmap color scale and data categories. Graphic 
    device must be initialized first.
}
\usage{
    bioMatrixLegend(heatmapNames=NULL, categoryNames=NULL, binaryNames=NULL, 
            heatmapMin=-3, heatmapMax=3, colorType="BlueWhiteRed")
}
\arguments{
    \item{heatmapNames}{
    character vector of length 2, name(s) of dataset for heatmap, e.g., 
    "RNASeq" and/or "miRNASeq"
    }
    \item{categoryNames}{
    character vector of length 2 or more, names of categories, e.g., 
    "Methylation High", "Methylation Low", ...
    }
    \item{binaryNames}{
    character vector of length 2 , names of binary data, e.g., 
    "DNA Amplification" and "DNA Deletion".
    }
    \item{heatmapMin}{
    numeric, minimum values of heatmap plot data, default -3 (z-scores)
    }
    \item{heatmapMax}{
    numeric, maximum values of heatmap plot data, default 3 (z-scores)
    }
    \item{colorType}{
    characte veector, one of "BlueWhiteRed", "GreenWhiteRed", "GreenYellowRed",
    "GreenBlackRed" , or "YellowToRed"
    }
}
\details{
    This function will plot legend on the bottom of matrix layout if any 
    argument is defined. The order of legend items (from left to right) is 
    heatmap color scale followed by colored boxes for category data legend 
    then colored points for binary data legend. 
}
\value{
    None 
}
\author{
    Henry Zhang
}
\examples{
    data(biomatrixPlotDemoData)
    plotBioMatrix(biomatrixPlotDemoData, summaryType="text")
    bioMatrixLegend(heatmapNames=c("RNASeq", "miRNASeq"), 
        categoryNames=c("Methyl H", "Methyl L"), 
        binaryNames=c("CN LOSS", "CN Gain"),   
        colorType="BlueWhiteRed")
}
\keyword{ methods }

