\name{methylDemoData}
\alias{methylDemoData}
\docType{data}
\title{
    Methylation Demo Data Set
}
\description{
    A data frame with beta values of 26 genes in 60 samples. Used for bar 
    plot on caOmicsV bioNetCircos layout demo, and for category plot on 
    caOmicsV bioNet layout demo after transforming to category data.
}
\usage{data("methylDemoData")}
\format{
    A data frame with 26 observations on the following 61 variables.
    \describe{
        \item{\code{Gene_Symbol}}{
            a factor with levels 
            \code{ACTN1} \code{ADAMTS13} \code{AMIGO3} \code{ATP2A1} 
            \code{BCO2} \code{CDKN3} \code{CFP} \code{CNBP} \code{COL15A1} 
            \code{CSRNP1} \code{CXCL12} \code{DBH} \code{DDX55} \code{ECM1} 
            \code{ELOVL1} \code{ESM1} \code{FAM81A} \code{FCN3} \code{KCNQ1} 
            \code{LEPREL1} \code{LIFR} \code{LILRA6} \code{LILRB5} 
            \code{LOC222699} \code{LOC283050} \code{LRRC16A} \code{LYVE1} 
            \code{MAEL} \code{MAN1B1} \code{MRPS25} \code{MT1F} \code{NIPA2} 
            \code{NPHP4} \code{NR5A2} \code{OLFML3} \code{PLVAP} \code{PROX1} 
            \code{PTH1R} \code{RBL2} \code{RCAN1} \code{RND3} \code{SEMA3F} 
            \code{SLC26A6} \code{TOMM40L} \code{TSEN34} \code{VDAC3}
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
    data(methylDemoData)
}
\keyword{datasets}
