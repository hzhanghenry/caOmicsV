\name{RNASeqDemoData}
\alias{RNASeqDemoData}
\docType{data}
\title{
    RNASeq Demo Data Set
}
\description{
    A data frame with RNASeq read counts of 26 genes in 60 samples. Used for 
    heatmap plot on caOimcsV bioNetCircos layout and bioMatirx layout after 
    transformed to log2 values
}
\usage{data("RNASeqDemoData")}
\format{
    A data frame with 26 observations on the following 61 variables.
    \describe{
        \item{\code{GeneSymbol}}{
            a factor with levels 
            \code{A1BG} \code{ADAMTS13} \code{AMIGO3} \code{ARMCX6} 
            \code{BCO2} \code{C17orf68} \code{CBLN4} \code{CDKN3} \code{CFP} 
            \code{COL13A1} \code{COL15A1} \code{CSRNP1} \code{CXCL12} 
            \code{DBH} \code{DNA2} \code{ECM1} \code{ESM1} \code{FAM20A} 
            \code{FCN3} \code{GJB2} \code{HOXD9} \code{KIF12} \code{LEPREL1} 
            \code{LIFR} \code{LILRB5} \code{LOC222699} \code{LOC441455} 
            \code{LYVE1} \code{MOBKL1A} \code{MT1F} \code{NPM3} \code{OLFML3} 
            \code{PCDHB7} \code{PLVAP} \code{PRICKLE4} \code{PTH1R} 
            \code{RCAN1} \code{RND3} \code{RNF125} \code{SEMA3F} \code{SIRT5} 
            \code{SLC26A6} \code{SPACA5} \code{TKTL1} \code{TOMM40L} 
            \code{UBE2J2} \code{ZNF273}
        }
        \item{\code{TCGA.BC.A216.Normal}}{a numeric vector}
        \item{\code{TCGA.BD.A2L6.Normal}}{a numeric vector}
        \item{\code{TCGA.BD.A3EP.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A113.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A114.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A118.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A119.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A11A.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A11B.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A11C.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A11D.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EB.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EC.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EG.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EH.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EI.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EJ.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EL.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A39V.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A39W.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A39X.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A39Z.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A3A1.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A3A2.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A3A3.Normal}}{a numeric vector}
        \item{\code{TCGA.EP.A12J.Normal}}{a numeric vector}
        \item{\code{TCGA.EP.A26S.Normal}}{a numeric vector}
        \item{\code{TCGA.ES.A2HT.Normal}}{a numeric vector}
        \item{\code{TCGA.FV.A23B.Normal}}{a numeric vector}
        \item{\code{TCGA.FV.A2QR.Normal}}{a numeric vector}
        \item{\code{TCGA.BC.A216.Tumor}}{a numeric vector}
        \item{\code{TCGA.BD.A2L6.Tumor}}{a numeric vector}
        \item{\code{TCGA.BD.A3EP.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A113.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A114.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A118.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A119.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A11A.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A11B.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A11C.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A11D.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EB.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EC.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EG.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EH.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EI.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EJ.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EL.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A39V.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A39W.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A39X.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A39Z.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A3A1.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A3A2.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A3A3.Tumor}}{a numeric vector}
        \item{\code{TCGA.EP.A12J.Tumor}}{a numeric vector}
        \item{\code{TCGA.EP.A26S.Tumor}}{a numeric vector}
        \item{\code{TCGA.ES.A2HT.Tumor}}{a numeric vector}
        \item{\code{TCGA.FV.A23B.Tumor}}{a numeric vector}
        \item{\code{TCGA.FV.A2QR.Tumor}}{a numeric vector}
    }
}
\value{
    None 
}
\examples{
    data(RNASeqDemoData)
}
\keyword{datasets}
