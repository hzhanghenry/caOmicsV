\name{getRelatedPlotData}
\alias{getRelatedPlotData}
\title{
    Extract subset from A Data Frame Based on Relational Information
}
\description{
    Extract subset of plot data based on relational information, e.g., 
    expression of miRNA or DNA copy number variations that are correlated to 
    differentially expressed genes 
}
\usage{
    getRelatedPlotData(omicsData, linkData, geneNames)
}
\arguments{
    \item{omicsData}{
    a data frame, the dataset from which subset is extracted
    }
    \item{linkData}{
    a data frame, usually gene names and their related items. The first column 
    must be the items to which the second item is linked to.
    }
    \item{geneNames}{
    character vector, subset of gene names for subset.
    }
}
\value{
    A data frame with subset of the first input data.
}
\author{
    Henry Zhang
}
\examples{
    data(miRNA)
    data(RNA2miRNA)
    geneNames <- as.character(RNA2miRNA[,1])
    miRNAexpr <- getRelatedPlotData(miRNA, RNA2miRNA, geneNames)
}
\keyword{ methods }
