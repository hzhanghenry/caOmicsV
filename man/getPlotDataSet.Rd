\name{getPlotDataSet}
\alias{getPlotDataSet}
\title{
    Prepare Data Set for caOmicsV Plot
}
\description{
    This function will validate each dataset then convert them to matrix and 
    wrap all of them in one list object.
}
\usage{
    getPlotDataSet(sampleNames, geneNames, sampleData, heatmapData=list(), 
        categoryData=list(), binaryData=list(), summaryData=list(), 
        secondGeneNames=NULL)
}
\arguments{
    \item{sampleNames}{
    character vector, sample names, must be same or exist in every data set.
    }
    \item{geneNames}{
    character vector, gene names, must be same or exist in every data set.
    }
    \item{sampleData}{
    data frame with rows for samples and columns for features.
    }
    \item{heatmapData}{
    list of data frame(s) for heatmap plot. The first column of each data 
    frame is row names and others are numeric values. The list could be empty, 
    or having one or more data frame in a list object. Heatmap data should be 
    log2 values or z-scores.
    }
    \item{categoryData}{
    list of data frame(s). The first column of each data frame is row names 
    and others are numeric values. The list could be empty, or having one or 
    more data frame in a list object
    }
    \item{binaryData}{
    list of data frame(s). The first column of each data frame is row names 
    and others are binary values. The listcould be none, or one or more data 
    frame in a list object
    }
    \item{summaryData}{
    list of data frames with summary information for samples (columns) or for 
    genes (rows). The first column is for ID following by one or more columns 
    of summary data.
    }
    \item{secondGeneNames}{
    character vector, gene names that will be plot on right side of biomatrix 
    plot layout
    }
}
\value{
    \item{sampleNames }{character verctor, sample names}
    \item{geneNames }{character verctor, gene names}
    \item{secondGeneNames }{character verctor, for example, miRNA names }
    \item{sampleInfo}{ 
    a data frame, sample information such as Tumor/Normal, age, diagnosis
    }
    \item{heatmapData}{
    list of data matrix(s), e.g., RNASeq read counts at gene level
    }
    \item{categoryData}{
    list of data matrix(s), such as SNP in a gene, homozygous, or 
    heterozygous, or wildtype
    }
    \item{binaryData}{
    list of data matrix(s), e.g., mutation status of the gene
    }
    \item{
    summaryInfo}{list of data matrix(s), such as percentage of highly 
    expressed miRNA in all samples
    }
}

\author{
    Henry Zhang
}
\examples{

    data(sampleDemoData) 
    data(RNA2miRNA)
    data(RNASeqDemoData)
    data(miRNADemoData)
    data(methylDemoData)
    data(CNVDemoData)

    sampleNames <- as.character(sampleDemoData[,1])
    geneNames <- as.character(RNA2miRNA[,1])
    secondGeneNames <- as.character(RNA2miRNA[,2])

    normals <- grep("Normal", colnames(RNASeqDemoData))
    tumors  <- grep("Tumor", colnames(RNASeqDemoData))
    tumorExpr <- RNASeqDemoData[, tumors]
    normalExpr <- RNASeqDemoData[, normals]
    meanLog2Fold <- log2(rowMeans(tumorExpr/normalExpr))
    summaryData <- data.frame(geneNames, meanLog2Fold)

    plotData <- getPlotDataSet(sampleNames, geneNames, sampleDemoData, 
        heatmapData=list(RNASeqDemoData, miRNADemoData), 
        categoryData=list(methylDemoData), 
        binaryData=list(CNVDemoData), 
        summaryData=list(summaryData), 
        secondGeneNames)
}
\keyword{ methods }

