\name{sortClinicalData}
\alias{sortClinicalData}
\title{
    Sort Clinical Data by a Column
}
\description{
    Sort sample information by a column defined by byItem (a column header). 
    This function is for sorting sample information only. 
}
\usage{
    sortClinicalData(clinicalData, byItem)
}
\arguments{
    \item{clinicalData}{
    A data frame with rows for samples and columns for features. Sample names 
    must be in the first column.
    }
    \item{byItem}{
    character vector of a feature (column header) by which the data will be 
    sorted.
    }
}
\value{
    \item{sampleData}{
    copy of the first argument with new row order.
    }
}
\author{
    Henry Zhang
}
\examples{
    data(sampleDemoData)
    sampleInfo <- sortClinicalData(sampleDemoData, colnames(sampleDemoData)[1])
}
\keyword{ methods }
