\name{getBioMatrixParameters}
\alias{
    getBioMatrixBasePositions
    getBioMatrixColumnPadding
    getBioMatrixDataAreaWidth
    getBioMatrixGeneLabelWidth
    getBioMatrixGeneNumber
    getBioMatrixLegendHeight
    getBioMatrixPhenotypeNumber
    getBioMatrixPlotAreaHeigth
    getBioMatrixPlotAreaWidth
    getBioMatrixRemarkWidth
    getBioMatrixSummaryWidth
    getBioMatrixRowPadding
    getBioMatrixSampleHeight
    getBioMatrixSampleIDHeight
    getBioMatrixSampleNumber
    getBioMatrixSampleWidth
}
\title{
    Methods to Get caOmicsV BioMatrix Plot Parameters
}
\description{
    Get methods to retrieve parameters for caOmicsV bioMatrix layout plot 
    stored in caOmicsV environment. bioMatrix layout must be initialized first.
}
\usage{
    getBioMatrixBasePositions()
    getBioMatrixColumnPadding()
    getBioMatrixDataAreaWidth()
    getBioMatrixGeneLabelWidth()
    getBioMatrixGeneNumber()
    getBioMatrixLegendHeight()
    getBioMatrixPhenotypeNumber()
    getBioMatrixPlotAreaHeigth()
    getBioMatrixPlotAreaWidth()
    getBioMatrixRemarkWidth()
    getBioMatrixSummaryWidth()
    getBioMatrixRowPadding()
    getBioMatrixSampleHeight()
    getBioMatrixSampleIDHeight()
    getBioMatrixSampleNumber()
    getBioMatrixSampleWidth()
}
\value{
    getBioMatrixBasePositions() returns a numeric matrix of defult x and y 
    coordinates of rectangles for each samples at a row.  

    getBioMatrixColumnPadding() returns a non-negative numeric value in inch 
    for padding between two samples, default 0.025.  

    getBioMatrixDataAreaWidth() returns a non-negative numeric value in inch 
    for width of data plot area, total samples times the sum of sample width 
    and columnPadding.  

    getBioMatrixDataRowTop() returns y coordinate for the top of a sample row.

    getBioMatrixGeneLabelWidth() returns a non-negative numeric value in inch 
    for length of left labels (gene names).  

    getBioMatrixGeneNumber() returns total number of genes to be plotted.  

    getBioMatrixLegendHeight() returns a non-negative numeric value in inch 
    for the height of legend area.

    getBioMatrixPhenotypeNumber() returns total number of phenotypes.  

    getBioMatrixPlotAreaHeigth() returns a non-negative numeric value in inch 
    representing the height of all plot areas (sample name area, data plot 
    area, and legend area).  

    getBioMatrixPlotAreaWidth() returns a non-negative numeric value in inch 
    for width of all plot areas (left labels (gene names), data plot area, and 
    right labels).  

    getBioMatrixRemarkWidth() returns a non-negative numeric value in inch for 
    width on the right side of data plot, usually for second gene labels.

    getBioMatrixSummaryWidth() returns a non-negative numeric value in inch  
    for width on the right side of data plot area to plot summary data plot.  

    getBioMatrixRowPadding() returns a non-negative numeric value in inch for 
    height of padding between two rows. 

    getBioMatrixSampleHeight() returns a non-negative numeric value in inch 
    for height of a sample row. 

    getBioMatrixSampleIDHeight() returns a non-negative numeric value in inch 
    for height of sample labels (on the top of phenotype plot area).  

    getBioMatrixSampleNumber() returns the total number of samples in plot 
    datasets.  

    getBioMatrixSampleWidth() returns a non-negative numeric value in inch for 
    width of a rectangle (sample).
}
\author{
    Henry Zhang
}
\examples{
    initializeBioMatrixPlot()

    positions      <- getBioMatrixBasePositions()
    colPadding     <- getBioMatrixColumnPadding()
    dataAreaWidth  <- getBioMatrixDataAreaWidth()
    geneNameWidth  <- getBioMatrixGeneLabelWidth()
    numOfGenes     <- getBioMatrixGeneNumber()
    legendHeight   <- getBioMatrixLegendHeight()
    numOfFeatures  <- getBioMatrixPhenotypeNumber()
    dataAreaHeight <- getBioMatrixPlotAreaHeigth()
    plotAreaWidth  <- getBioMatrixPlotAreaWidth()
    sumAreaWidth   <- getBioMatrixRemarkWidth()
    rowPadding     <- getBioMatrixRowPadding()
    sampleHeight   <- getBioMatrixSampleHeight()
    sampleIDHeight <- getBioMatrixSampleIDHeight()
    numOfSamples   <- getBioMatrixSampleNumber()
    sampleWidth    <- getBioMatrixSampleWidth()
}
\keyword{ methods }

