\name{sortOmicsDataByRow}
\alias{sortOmicsDataByRow}
\title{
    Sort Omics Data by Row
}
\description{
    Sort omics data by row (genes) in order to get omics data in a specific 
    order, e.g., ordered by p values or by fold changes in expression data. 
    After sorting, the order of row names in omics data will be same as the 
    second argument.
}
\usage{
    sortOmicsDataByRow(omicsData, geneNames)
}
\arguments{
    \item{omicsData}{
    a data frame that holds genomic data such as gene expression, SNV, RNASeq. 
    The first column must be gene names that are same as the second argument 
    except of order.
    }
    \item{geneNames}{
    character vector, gene names in a given order (such as by p values). 
    Redundant gene names are allowed.
    }
}
\value{
    \item{omicsData}{
    copy of the first argument with new row order.
    }
}
\author{
    Henry Zhang
}
\examples{
    data(RNASeqDemoData)
    geneNames <- as.character(RNASeqDemoData[,1])
    geneNames <- geneNames[order(RNASeqDemoData[,2])]
    expr <- sortOmicsDataByRow(RNASeqDemoData, geneNames)
}
\keyword{ methods }
