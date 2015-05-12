\name{sortOmicsDataByColumn}
\alias{sortOmicsDataByColumn}
\title{
    Sort Omics Data by Column Header
}
\description{
    Sort omics data by column header (sample names) based on the order of the 
    second argument. The sample names in the first argument must have same 
    items as the second argument except of order. After sorting, the sample 
    names in omics data and the second argument will be in same order.
}
\usage{
sortOmicsDataByColumn(omicsData, sampleNames)
}
\arguments{
    \item{omicsData}{
    A data frame that holds genomic data such as gene expression, SNV, RNASeq 
    ... The column headers must be the sample names that are same as the 
    sample names in clinical data.
    }
    \item{sampleNames}{
    character vector, sample names in a given order (such as diagnosis). No 
    redundant entries allowed.
    }
}
\value{
    \item{omicsData}{
    copy of the first argument with columns in new order.
    }
}
\author{
    Henry Zhang
}
\examples{
    data(RNASeqDemoData)
    sampleNames <- colnames(RNASeqDemoData)[-1]
    sampleNames <- sampleNames[length(sampleNames):1]
    expr <- sortOmicsDataByColumn(RNASeqDemoData, sampleNames)
}
\keyword{ methods }
