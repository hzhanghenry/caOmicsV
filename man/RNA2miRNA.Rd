\name{RNA2miRNA}
\alias{RNA2miRNA}
\docType{data}
\title{
    Paired miRNA and Genes with Negative Correlation Coefficients
}
\description{
    A data frame contains 26 genes and their most negative correlated miRNA. 
    Used for demo of generation a new data frame with gene names and related 
    miRNA read counts.
}
\usage{data("RNA2miRNA")}
    \format{
        A data frame with 26 observations on the following 2 variables.
        \describe{
        \item{\code{GeneSymbol}}{
        a factor with levels 
        \code{ADAMTS13} \code{AMIGO3} \code{BCO2} \code{CDKN3} \code{CFP} 
        \code{COL15A1} \code{CSRNP1} \code{CXCL12} \code{DBH} \code{ECM1} 
        \code{ESM1} \code{FCN3} \code{LEPREL1} \code{LIFR} \code{LILRB5} 
        \code{LOC222699} \code{LYVE1} \code{MT1F} \code{OLFML3} \code{PLVAP} 
        \code{PTH1R} \code{RCAN1} \code{RND3} \code{SEMA3F} \code{SLC26A6} 
        \code{TOMM40L}
        }
        \item{\code{miRNA_ID}}{
        a factor with levels 
        \code{hsa.mir.101.1} \code{hsa.mir.101.2} \code{hsa.mir.10b} 
        \code{hsa.mir.139} \code{hsa.mir.424} \code{hsa.mir.450b}
        }
    }
}
\value{
    None 
}
\examples{
    data(RNA2miRNA)
}
\keyword{datasets}
