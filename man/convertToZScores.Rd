\name{convertToZScores}
\alias{convertToZScores}
\title{
    Calculate z-scores for A Data Matrix
}
\description{
    Calculate z-scores for data values in a data frame. The row ID must be in 
    the first column followed by data values.
}
\usage{
    convertToZScores(exprData)
}
\arguments{
    \item{exprData}{
    A data frame with first column as row IDs and others are numeric values.
    }
}
\value{
    A data frame with z scores for each row. The first column is still row IDs.
}
\author{
    Henry Zhang
}
\examples{
    data(RNASeqDemoData)
    exprZ <- convertToZScores(RNASeqDemoData)
}
\keyword{ methods }

