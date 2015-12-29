\name{plotBioMatrix}
\alias{plotBioMatrix}
\title{
    Default Plot Method for caOmicsV bioMatrix Layout Plot
}
\description{
    A sample way to plot sample information and genomic data with default 
    settings.
}
\usage{
    plotBioMatrix(dataSet, summaryType=c("text", "bar"), summarybyRow=TRUE, 
        heatmapMax=NULL, heatmapMin=NULL, heatmapColor="BlueWhiteRed")
}
\arguments{
    \item{dataSet}{
    an object, returned from method of getPlotDataSet(), with all plot data
    }
    \item{summaryType}{
    character vector, either "text" or "bar",  for plot type of summary data
    }
    \item{summarybyRow}{
    logic, if the summary data is for each row or for each column
    }
    \item{heatmapMax}{
    numeric, maximum value for heatmap plot, set to NULL to use the maximum 
    value in input data.
    }
    \item{heatmapMin}{
    numeric, minimum value for heatmap plot, set to NULL to use the minimum 
    value in input data.
    }
    \item{heatmapColor}{
    character vector,one of "BlueWhiteRed", "GreenWhiteRed", "GreenYellowRed", 
    "GreenBlackRed" , or "YellowToRed"
    }
}
\value{
    None
}
\author{
    Henry Zhang
}
\examples{
    data(biomatrixPlotDemoData)
    plotBioMatrix(biomatrixPlotDemoData, summaryType="text", 
        summarybyRow=TRUE, heatmapMax=3, heatmapMin=-3, 
        heatmapColor="BlueWhiteRed") 
}
\keyword{ methods }
