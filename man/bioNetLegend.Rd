\name{bioNetLegend}
\alias{bioNetLegend}
\title{
    Draw Legend for caOmicsV bioNet Plot
}
\description{
    The bioNet legend includes a heatmap color scale and names for each track, 
    such as "1. Tissue: T(red), N(blue)", "2. Methylation", "3. miRNA 
    hsa-mir-424", "4. Gene Expression", ...  A graphic device and igraph 
    object must be initialized first.
}
\usage{
    bioNetLegend(dataNames, textCoor=NULL, heatmapCoor=NULL,
        scaleWidth, scaleHeight, heatmapMin=-3, heatmapMax=3,  
        colorType="BlueWhiteRed", direction="h")
}
\arguments{
    \item{dataNames}{
    character vector, names of data on each circular track
    }
    \item{textCoor}{
    numeric vector of length 2, x and y coordinates for legend text. 
    }
    \item{heatmapCoor}{
    numeric vector of length 2, x and y coordinates for heatmap colour scale
    }
    \item{scaleWidth}{
    non-negative numeric, length (width) of heatmap color scale
    }
    \item{scaleHeight}{
    non-negative numeric, height of heatmap color scale
    }
    \item{heatmapMin}{
    numeric, minimum value of heatmap color scale, default -3 (z-score)
    }
    \item{heatmapMax}{
    numeric, maximum value of heatmap color scale, default 3 (z-score)
    }
    \item{colorType}{
    character vector, one of "BlueWhiteRed", "GreenWhiteRed", "GreenYellowRed", 
    "GreenBlackRed" , or "YellowToRed"
    }
    \item{direction}{
    character, direction of heatmap color scale, either 'h' for horizontal or 
    'v' for vertical.
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
    dataNames <- c("Tissue Type", "RNASeq", "miRNASeq", "Methylation", "CNV")
    bioNetLegend(dataNames)
}
\keyword{ methods }
