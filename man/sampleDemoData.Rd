\name{sampleDemoData}
\alias{sampleDemoData}
\docType{data}
\title{
    Sample Tissue Type Demo Data Set
}
\description{
    A data frame with sample names and tissue types for 60 samples. Used to 
    plot tissue type on both caOmicsV bioMatrix and bioNetCircos layout.
}
\usage{data("sampleDemoData")}
\format{
    A data frame with 60 observations on the following 2 variables.
    \describe{
        \item{\code{sample_names}}{
            a factor with levels 
                \code{TCGA.BC.A216.Normal} \code{TCGA.BC.A216.Tumor} 
                \code{TCGA.BD.A2L6.Normal} \code{TCGA.BD.A2L6.Tumor} 
                \code{TCGA.BD.A3EP.Normal} \code{TCGA.BD.A3EP.Tumor} 
                \code{TCGA.DD.A113.Normal} \code{TCGA.DD.A113.Tumor} 
                \code{TCGA.DD.A114.Normal} \code{TCGA.DD.A114.Tumor} 
                \code{TCGA.DD.A118.Normal} \code{TCGA.DD.A118.Tumor} 
                \code{TCGA.DD.A119.Normal} \code{TCGA.DD.A119.Tumor} 
                \code{TCGA.DD.A11A.Normal} \code{TCGA.DD.A11A.Tumor} 
                \code{TCGA.DD.A11B.Normal} \code{TCGA.DD.A11B.Tumor} 
                \code{TCGA.DD.A11C.Normal} \code{TCGA.DD.A11C.Tumor} 
                \code{TCGA.DD.A11D.Normal} \code{TCGA.DD.A11D.Tumor} 
                \code{TCGA.DD.A1EB.Normal} \code{TCGA.DD.A1EB.Tumor} 
                \code{TCGA.DD.A1EC.Normal} \code{TCGA.DD.A1EC.Tumor} 
                \code{TCGA.DD.A1EG.Normal} \code{TCGA.DD.A1EG.Tumor} 
                \code{TCGA.DD.A1EH.Normal} \code{TCGA.DD.A1EH.Tumor} 
                \code{TCGA.DD.A1EI.Normal} \code{TCGA.DD.A1EI.Tumor} 
                \code{TCGA.DD.A1EJ.Normal} \code{TCGA.DD.A1EJ.Tumor} 
                \code{TCGA.DD.A1EL.Normal} \code{TCGA.DD.A1EL.Tumor} 
                \code{TCGA.DD.A39V.Normal} \code{TCGA.DD.A39V.Tumor} 
                \code{TCGA.DD.A39W.Normal} \code{TCGA.DD.A39W.Tumor} 
                \code{TCGA.DD.A39X.Normal} \code{TCGA.DD.A39X.Tumor} 
                \code{TCGA.DD.A39Z.Normal} \code{TCGA.DD.A39Z.Tumor} 
                \code{TCGA.DD.A3A1.Normal} \code{TCGA.DD.A3A1.Tumor} 
                \code{TCGA.DD.A3A2.Normal} \code{TCGA.DD.A3A2.Tumor} 
                \code{TCGA.DD.A3A3.Normal} \code{TCGA.DD.A3A3.Tumor} 
                \code{TCGA.EP.A12J.Normal} \code{TCGA.EP.A12J.Tumor} 
                \code{TCGA.EP.A26S.Normal} \code{TCGA.EP.A26S.Tumor} 
                \code{TCGA.ES.A2HT.Normal} \code{TCGA.ES.A2HT.Tumor} 
                \code{TCGA.FV.A23B.Normal} \code{TCGA.FV.A23B.Tumor} 
                \code{TCGA.FV.A2QR.Normal} \code{TCGA.FV.A2QR.Tumor}
        }
        \item{\code{tissue_type}}{
                a factor with levels
                \code{Normal} \code{Tumor}
        }
    }
}
\value{
    None 
}
\examples{
    data(sampleDemoData)
}
\keyword{datasets}
