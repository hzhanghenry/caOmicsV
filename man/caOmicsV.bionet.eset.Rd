\name{caOmicsV.bionet.eset}
\alias{caOmicsV.bionet.eset}
\docType{data}
\title{
    Demo Dataset for caOmicsV bioNetCircos Plot
}
\description{
    eset for demo of caOmicsV bionetCircos Plot
}
\usage{data("caOmicsV.bionet.eset")}
\format{
    The format is: List of 8 

    $ sampleNames: chr [1:60] "BC.A216.Normal" "BD.A2L6.Normal" 
            "BD.A3EP.Normal" "DD.A113.Normal" ...

    $ geneNames: chr [1:26] "ECM1" "SLC26A6" "ADAMTS13" "FCN3" ...

    $ secondGeneNames: NULL 

    $ sampleInfo: chr [1:2, 1:60] "TCGA.BC.A216.Normal"   
        "Solid Tissue Normal" "TCGA.BD.A2L6.Normal" "Solid Tissue Normal" ...
        ..- attr(*, "dimnames")=List of 2
        .. ..$ : chr [1:2] "sampleID" "sample_type"
        .. ..$ : chr [1:60] "TCGA.BC.A216.Normal" "TCGA.BD.A2L6.Normal" 
            "TCGA.BD.A3EP.Normal" "TCGA.DD.A113.Normal" ...

    $ heatmapData:List of 2
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
        ..$ : num [1:26, 1:60] 0.825 0.364 0.545 0.798 0.606 ...
        .. ..- attr(*, "dimnames")=List of 2
        .. .. ..$ : chr [1:26] "ECM1" "SLC26A6" "ADAMTS13" "FCN3" ...
        .. .. ..$ : chr [1:60] "TCGA.BC.A216.Normal" "TCGA.BD.A2L6.Normal" 
            "TCGA.BD.A3EP.Normal" "TCGA.DD.A113.Normal" ...

    $ binaryData     :List of 1 
        ..$ : num [1:26, 1:60] 0.0266 0.0069 0.0193 0.0117 0.0081 0.0242 
            0.0069 0.0079 0.0348 0.01 ...
        .. ..- attr(*, "dimnames")=List of 2
        .. .. ..$ : chr [1:26] "ECM1" "SLC26A6" "ADAMTS13" "FCN3" ...
        .. .. ..$ : chr [1:60] "TCGA.BC.A216.Normal" "TCGA.BD.A2L6.Normal" 
            "TCGA.BD.A3EP.Normal" "TCGA.DD.A113.Normal" ... 

    $ summaryInfo    : NULL
}
\value{
    None
}
\examples{
    data("caOmicsV.bionet.eset")
}
\keyword{datasets}
