\name{miRNA}
\alias{miRNA}
\docType{data}
\title{
    Selected miRNA Read Counts
}
\description{
    A data frame with read counts of 14 miRNA in 86 samples.
}
\usage{data("miRNA")}
\format{
    A data frame with 14 observations on the following 86 variables.
    \describe{
        \item{\code{miRNA_ID}}{
            a factor with levels
            \code{hsa.mir.101.1} \code{hsa.mir.101.2} \code{hsa.mir.10b} 
            \code{hsa.mir.1180} \code{hsa.mir.125b.2} \code{hsa.mir.139} 
            \code{hsa.mir.142} \code{hsa.mir.151} \code{hsa.mir.183} 
            \code{hsa.mir.22} \code{hsa.mir.25} \code{hsa.mir.424} 
            \code{hsa.mir.450b} \code{hsa.mir.93}
        }
        \item{\code{TCGA.2Y.A9H1.Tumor}}{a numeric vector}
        \item{\code{TCGA.BC.4073.Tumor}}{a numeric vector}
        \item{\code{TCGA.BC.A10W.Tumor}}{a numeric vector}
        \item{\code{TCGA.BC.A112.Tumor}}{a numeric vector}
        \item{\code{TCGA.BC.A216.Tumor}}{a numeric vector}
        \item{\code{TCGA.BC.A216.Normal}}{a numeric vector}
        \item{\code{TCGA.BD.A2L6.Tumor}}{a numeric vector}
        \item{\code{TCGA.BD.A2L6.Normal}}{a numeric vector}
        \item{\code{TCGA.BD.A3EP.Tumor}}{a numeric vector}
        \item{\code{TCGA.BD.A3EP.Normal}}{a numeric vector}
        \item{\code{TCGA.BW.A5NP.Tumor}}{a numeric vector}
        \item{\code{TCGA.CC.A1HT.Tumor}}{a numeric vector}
        \item{\code{TCGA.CC.A7IG.Tumor}}{a numeric vector}
        \item{\code{TCGA.CC.A9FS.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A113.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A113.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A114.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A114.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A118.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A118.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A119.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A119.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A11A.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A11A.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A11B.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A11B.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A11C.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A11C.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A11D.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A11D.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EB.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EB.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EC.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EC.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EE.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EG.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EG.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EH.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EH.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EI.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EI.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EJ.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EJ.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A1EL.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A1EL.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A39V.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A39V.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A39W.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A39W.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A39X.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A39X.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A39Z.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A39Z.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A3A1.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A3A1.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A3A2.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A3A2.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A3A3.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A3A3.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A3A6.Normal}}{a numeric vector}
        \item{\code{TCGA.DD.A4NG.Tumor}}{a numeric vector}
        \item{\code{TCGA.DD.A4NV.Tumor}}{a numeric vector}
        \item{\code{TCGA.ED.A4XI.Tumor}}{a numeric vector}
        \item{\code{TCGA.ED.A82E.Tumor}}{a numeric vector}
        \item{\code{TCGA.EP.A12J.Tumor}}{a numeric vector}
        \item{\code{TCGA.EP.A12J.Normal}}{a numeric vector}
        \item{\code{TCGA.EP.A26S.Tumor}}{a numeric vector}
        \item{\code{TCGA.EP.A26S.Normal}}{a numeric vector}
        \item{\code{TCGA.ES.A2HS.Tumor}}{a numeric vector}
        \item{\code{TCGA.ES.A2HT.Tumor}}{a numeric vector}
        \item{\code{TCGA.ES.A2HT.Normal}}{a numeric vector}
        \item{\code{TCGA.FV.A23B.Tumor}}{a numeric vector}
        \item{\code{TCGA.FV.A23B.Normal}}{a numeric vector}
        \item{\code{TCGA.FV.A2QR.Tumor}}{a numeric vector}
        \item{\code{TCGA.FV.A2QR.Normal}}{a numeric vector}
        \item{\code{TCGA.FV.A4ZQ.Tumor}}{a numeric vector}
        \item{\code{TCGA.G3.A3CH.Tumor}}{a numeric vector}
        \item{\code{TCGA.G3.A6UC.Tumor}}{a numeric vector}
        \item{\code{TCGA.G3.AAV3.Tumor}}{a numeric vector}
        \item{\code{TCGA.K7.A5RF.Tumor}}{a numeric vector}
        \item{\code{TCGA.LG.A9QD.Tumor}}{a numeric vector}
        \item{\code{TCGA.O8.A75V.Tumor}}{a numeric vector}
        \item{\code{TCGA.RC.A7SH.Tumor}}{a numeric vector}
        \item{\code{TCGA.UB.AA0U.Tumor}}{a numeric vector}
        \item{\code{TCGA.ZP.A9CZ.Tumor}}{a numeric vector}
    }
}
\value{
    None 
}
\examples{
    data(miRNA)
}
\keyword{datasets}
