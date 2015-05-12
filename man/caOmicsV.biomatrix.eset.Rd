\name{caOmicsV.biomatrix.eset}
\alias{caOmicsV.biomatrix.eset}
\docType{data}
\title{
    Demo Dataset for caOmicsV bioMatrix Plot
}
\description{
    eset for demo of caOmicsV bioMatrix Plot
}
\usage{data("caOmicsV.biomatrix.eset")}
\format{
    The format is: List of 8

    $ sampleNames: chr [1:60] "BC.A216.Normal" "BD.A2L6.Normal" 
            "BD.A3EP.Normal" "DD.A113.Normal" ...

    $ geneNames: chr [1:26] "ECM1" "SLC26A6" "ADAMTS13" "FCN3" ... 

    $ secondGeneNames: chr [1:26] "hsa-mir-10b" "hsa-mir-139" 
            "hsa-mir-10b" "hsa-mir-10b" ...

    $ sampleInfo: chr [1:2, 1:60] "TCGA.BC.A216.Normal"   
        "Solid Tissue Normal" "TCGA.BD.A2L6.Normal" "Solid Tissue Normal" ...
        ..- attr(*, "dimnames")=List of 2
        .. ..$ : chr [1:2] "sampleID" "sample_type"
        .. ..$ : chr [1:60] "TCGA.BC.A216.Normal" "TCGA.BD.A2L6.Normal" 
                "TCGA.BD.A3EP.Normal" "TCGA.DD.A113.Normal" ...

    $ heatmapData    :List of 2
        ..$ : num [1:26, 1:60] 1.157 -0.623 0.667 0.976 0.868 ...
        .. ..- attr(*, "dimnames")=List of 2
        .. .. ..$ : chr [1:26] "ECM1" "SLC26A6" "ADAMTS13" "FCN3" ...
        .. .. ..$ : chr [1:60] "TCGA.BC.A216.Normal" "TCGA.BD.A2L6.Normal" 
                "TCGA.BD.A3EP.Normal" "TCGA.DD.A113.Normal" ...
                
        ..$ : num [1:26, 1:60] -1.5 1.08 -1.5 -1.5 -1.5 ...
        .. ..- attr(*, "dimnames")=List of 2
        .. .. ..$ : chr [1:26] "ECM1" "SLC26A6" "ADAMTS13" "FCN3" ...
        .. .. ..$ : chr [1:60] "TCGA.BC.A216.Normal" "TCGA.BD.A2L6.Normal" 
                "TCGA.BD.A3EP.Normal" "TCGA.DD.A113.Normal" ...

    $ categoryData   :List of 1
        ..$ : num [1:26, 1:60] 1 0 1 1 1 0 0 0 1 0 ...
        .. ..- attr(*, "dimnames")=List of 2
        .. .. ..$ : chr [1:26] "ECM1" "SLC26A6" "ADAMTS13" "FCN3" ...
        .. .. ..$ : chr [1:60] "TCGA.BC.A216.Normal" "TCGA.BD.A2L6.Normal" 
                "TCGA.BD.A3EP.Normal" "TCGA.DD.A113.Normal" ...

    $ binaryData     :List of 2 
        ..$ : num [1:26, 1:60] 0 0 0 0 0 0 0 0 0 0 ... 
        .. ..- attr(*, "dimnames")=List of 2 
        .. .. ..$ : chr [1:26] "ECM1" "SLC26A6" "ADAMTS13" "FCN3" ... 
        .. .. ..$ : chr [1:60] "TCGA.BC.A216.Normal" "TCGA.BD.A2L6.Normal" 
                "TCGA.BD.A3EP.Normal" "TCGA.DD.A113.Normal" ... 

        ..$ : num [1:26, 1:60] 0 0 0 0 0 0 0 0 0 0 ...
        .. ..- attr(*, "dimnames")=List of 2
        .. .. ..$ : chr [1:26] "ECM1" "SLC26A6" "ADAMTS13" "FCN3" ...
        .. .. ..$ : chr [1:60] "TCGA.BC.A216.Normal" "TCGA.BD.A2L6.Normal" 
                "TCGA.BD.A3EP.Normal" "TCGA.DD.A113.Normal" ...

    $ summaryInfo    :List of 1
        ..$ : chr [1:26, 1:2] "ECM1" "SLC26A6" "ADAMTS13" "FCN3" ...
        .. ..- attr(*, "dimnames")=List of 2
        .. .. ..$ : chr [1:26] "ECM1" "SLC26A6" "ADAMTS13" "FCN3" ...
        .. .. ..$ : chr [1:2] "GeneSymbol" "logFC"
}
\value{
    None
}
\examples{
    data("caOmicsV.biomatrix.eset")
}
\keyword{datasets}
