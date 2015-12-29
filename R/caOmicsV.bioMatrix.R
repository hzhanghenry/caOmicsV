#
#    caOmicsV.bioMatrixPlot
#
#    Prerequisites: 
#
#    All datasets or arguments passed to functions have to have same samples, 
#    i.e., number of samples and order of samples must be same
#
#    Data types supported and desired graphic patterns:
#
#    Continuous data:
#
#        Gene/RNA expression: log2, plotted as heatmap or colored rectangle
#        miRNA expression:    log2, plotted as heatmap or colored rectangle
#
#     Categorical data:
#
#        Phenotypes: tumor, normal, ...,  plotted as colored rectangle  
#        CNVs: deletion, normal, amplification, plotted as colored points
#        Methylation: high, low, none, plotted as colored outline
#
#    Binary data:
#
#        Mutations/SNPs: plotted as colored points, positive only 
#        DNA indels:     plotted as colored points, positive only 
#        Gene fusions:   plotted as colored points, positive only 
#
#    Useful graphic parameters:
#
#    cex    number indicating the amount by which plotting text and symbols   
#            should be scaledrelative to the default. 1=default, 1.5 is 50% 
#            larger, 0.5 is 50% smaller, etc.
# 
#    ps    font point size (roughly 1/72 inch),  text size will be ps*cex
#
#    __________________________________________________________________________
#    **************************************************************************






#    __________________________________________________________________________
#    **************************************************************************
#
#    Initialize matrix plot. Call this function when new layout is needed. 
#
#    Arguments:
#
#    numOfGenes:       non-negative integer, number of genes
#    numOfSamples:     non-negative integer, number of samples
#    numOfPhenotypes:  non-negative integer, number of phenotypes
#
#    sampleHeight:     non-negative numeric, height of rectangle for each 
#                      sample in inch, default is 0.4
#    sampleWidth:      non-negative numeric, width of rectangle for each 
#                      sample in inch, default is 0.1
#    columnPadding:    non-negative numeric, padding between two samples 
#                      in inch, default is 0.025
#    rowPadding:       non-negative numeric, padding between two rows in 
#                      inch, default is 0.1
#    geneNameWith:    non-negative numeric, width of gene labelling in inch, 
#                      default is 1
#    remarkWidth:     non-negative numeric, width of remark area in inch, 
#                      default is 1
#    sampleNameHeight: non-negative numeric, height of column label area
#    legendHeight:     non-negative numeric, height of legend in inch   
#
#    Returned value:    None
#
#    Example: caOmicsV.InitializeBioMatrixPlot(numOfGenes=100,  
#                    numOfSamples=100, numOfPhenotypes=1, sampleHeight=0.4, 
#                    sampleWidth=0.1, columnPadding=0.025, rowPadding=0.1, 
#                    geneNameWidth=1, remarkWidth=1, 
#                    sampleNameHeight=1,legendHeight=1)
#
#    Last revised on September 3, 2014
#

initializeBioMatrixPlot <- function(numOfGenes=100, numOfSamples=100,
        numOfPhenotypes=1, sampleHeight=0.4, sampleWidth=0.1, 
        columnPadding=0.025, rowPadding=0.1, geneNameWidth=1,
        remarkWidth=1, summaryWidth=1, sampleNameHeight=1, legendHeight=1) {

    if(is.numeric(numOfGenes) == FALSE || 
        is.numeric(numOfSamples) == FALSE || 
        is.numeric(numOfPhenotypes) == FALSE || 
        is.numeric(sampleHeight) == FALSE || 
        is.numeric(sampleWidth) == FALSE || 
        is.numeric(columnPadding) == FALSE ||
        is.numeric(rowPadding) == FALSE  || 
        is.numeric(geneNameWidth) == FALSE || 
        is.numeric(remarkWidth) == FALSE || 
        is.numeric(summaryWidth) == FALSE ||
        is.numeric(sampleNameHeight) == FALSE || 
        is.numeric(legendHeight) == FALSE ) {
            stop("All arguments must be non-negative numeric.\n") 
    }

    if(numOfGenes<0 || numOfSamples<0 || numOfPhenotypes<0 || 
        sampleHeight<0 || sampleWidth<0 || columnPadding<0 || 
        rowPadding<0 || geneNameWidth<0 || remarkWidth<0 || 
        sampleNameHeight<0 || legendHeight<0 || summaryWidth<0)  {
            stop("Negative value is not allowed for argument(s)") 
    }

    setBioMatrixPlotParameters(numOfGenes, numOfSamples, 
            numOfPhenotypes, sampleHeight, sampleWidth, columnPadding, 
            rowPadding, geneNameWidth, remarkWidth, sampleNameHeight, 
            summaryWidth, legendHeight)
    setBioMatrixBaseCoordinates(numOfSamples, sampleWidth, columnPadding, 
            sampleHeight, geneNameWidth)
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Put biomatrix plot parameters to CA_OMICS_ENV environment. This function
#    is for internal use only and all arguments are validated outside first.
#
#    Arguments:
#
#        numOfGenes:      non-negative integer, number of genes to plot
#        numOfSamples:    non-negative integer, number of samples to plot
#        numOfPhenotypes: non-negative integer, number of phenotypes to plot
#
#        sampleHeight:     non-negative numeric, height of rectangle for each
#                          sample in inch, default 0.4
#        sampleWidth:      non-negative numeric, width of rectangle for each
#                          sample in inch, default 0.1
#        columnPadding:    non-negative numeric, padding between two samples
#                          in inch, default 0.025
#        rowPadding:       non-negative numeric, padding between two samples
#                          in inch, default 0.1
#        geneNameWith:     non-negative numeric, width of gene labelling in
#                          inch, default 1
#        remarkWidth:      non-negative numeric, width of remark area in inch
#                          default 1
#        sampleNameHeight: non-negative numeric, height of rectangle for
#                          sample names in inch, default 1
#        legendHeight:      non-negative numeric, height of legend in inch
#
#    Example:   Internal use.
#    Last updated on September 4, 2014
#

setBioMatrixPlotParameters <- function(numOfGenes, numOfSamples, 
    numOfPhenotypes, sampleHeight, sampleWidth, columnPadding, rowPadding,
    geneNameWidth, remarkWidth, summaryWidth, sampleNameHeight, legendHeight) {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    caOmicsVEnvironment[["BioMatrix_Total_Genes"]]      <- numOfGenes
    caOmicsVEnvironment[["BioMatrix_Total_Samples"]]    <- numOfSamples
    caOmicsVEnvironment[["BioMatrix_Total_Phenotypes"]] <- numOfPhenotypes

    caOmicsVEnvironment[["BioMatrix_Sample_Height"]]    <- sampleHeight
    caOmicsVEnvironment[["BioMatrix_Sample_Width"]]     <- sampleWidth
    caOmicsVEnvironment[["BioMatrix_Column_Padding"]]   <- columnPadding
    caOmicsVEnvironment[["BioMatrix_Row_Padding"]]      <- rowPadding

    caOmicsVEnvironment[["BioMatrix_Label_Width"]]      <- geneNameWidth
    caOmicsVEnvironment[["BioMatrix_Remark_Width"]]     <- remarkWidth
    caOmicsVEnvironment[["BioMatrix_Summary_Width"]]    <- summaryWidth
    caOmicsVEnvironment[["BioMatrix_SampleID_Height"]]  <- sampleNameHeight
    caOmicsVEnvironment[["BioMatrix_Legend_Height"]]    <- legendHeight

    dataAreaWidth  <- numOfSamples*(sampleWidth+columnPadding)
    caOmicsVEnvironment[["BioMatrix_DataArea_Width"]]   <- dataAreaWidth

    setCaOmicsVColors();
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    Set up x and y coordinates for one row of phenotype or omics data value  
#    plot. The x-coordinates will have no change at any time and y coordinates 
#    will be modified based on the row number before plot. This function is 
#    for internal use only and all arguments are validated outside in advance.
#
#    Arguments:
#
#        numOfSamples:  non-negative integer, number of samples to be plotted
#        sampleWidth:   non-negative numeric, width of rectangle for each 
#                       sample in inch, default 0.1
#        columnPadding: non-negative numeric, padding between two samples in 
#                       inch, default 0.025
#        sampleHeight:  non-negative numeric, height of rectangle for each  
#                       sample in inch, default 0.4
#        geneNameWith:  non-negative numeric, width of gene labelling in inch, 
#                       default 1
#
#   Example:    Internal use
#   Last updated on September 4, 2014
#

setBioMatrixBaseCoordinates <- function(numOfSamples, sampleWidth, 
            columnPadding, sampleHeight, geneNameWidth) {

    lastLeft  <- geneNameWidth + (numOfSamples-1)*(sampleWidth+columnPadding)
    lastRight <- lastLeft + sampleWidth
    sampleAreaWidth <- sampleWidth+columnPadding

    xleft   <- seq(from=geneNameWidth, to=lastLeft, by=sampleAreaWidth)
    xright  <- seq(from=geneNameWidth+sampleWidth, to=lastRight, 
                        by=sampleAreaWidth)
    ytop    <- rep(0, numOfSamples)
    ybottom <- rep(sampleHeight, numOfSamples)

    basePositions <- cbind(xleft, ybottom, xright, ytop)
    colnames(basePositions) <- c("xleft", "ybottom", "xright", "ytop")

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())
    caOmicsVEnvironment[["BioMatrix_Base_Positions"]] <- basePositions
}






#    __________________________________________________________________________
#    **************************************************************************
#
#    Set up plot area including of sample name area (sampleHeight), phenotype 
#    area, gene label area (geneNameWidth), and remark area (remarkWidth, for 
#    legend and other descriptions). The sampleHeight,  geneNameWidth, and 
#    remarkWidth could be increased by user to make more space for additional
#    item plot
#
#    Prerequisite: InitializeBioMatrixPlot(...) must be called first
#
#    Argument:  None
#    Returned value:    None
#
#    Example:
#
#        InitializeBioMatrixPlot(numOfGenes=100, numOfSamples=100, 
#                numOfPhenotypes=1, sampleHeight=0.4, sampleWidth=0.1, 
#                columnPadding=0.025, rowPadding=0.1, geneNameWidth=1, 
#                remarkWidth=1, sampleNameHeight=1)
#        setBioMatrixPlotArea()
#
#    Last updated on: September 5, 2014
#
#

setBioMatrixPlotArea <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    totalGenes     <- getBioMatrixGeneNumber()
    totalPhenotypes <- getBioMatrixPhenotypeNumber()

    geneNameWith   <- getBioMatrixGeneLabelWidth()
    dataAreaWidth  <- getBioMatrixDataAreaWidth()
    remarkWidth    <- getBioMatrixRemarkWidth()
    summaryWidth   <- getBioMatrixSummaryWidth()

    sampleHeight   <- getBioMatrixSampleHeight()
    rowPadding     <- getBioMatrixRowPadding()
    nameHeight     <- getBioMatrixSampleIDHeight()
    legendHeight   <- getBioMatrixLegendHeight()

    plotWidth  <- geneNameWith + dataAreaWidth + remarkWidth + summaryWidth
    rowHeight  <- sampleHeight + rowPadding
    plotHeight <- (totalGenes+totalPhenotypes)*rowHeight + 
                            nameHeight + legendHeight

    plot.new()
    plot.window(xlim=c(0, plotWidth), ylim=c(0, plotHeight))
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Methods to retrieving of bioMatrixPlot objects in CA_OMICS_ENV environment
#

getBioMatrixGeneNumber <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_Total_Genes"]])
}

getBioMatrixSampleNumber <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_Total_Samples"]])
}

getBioMatrixPhenotypeNumber <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_Total_Phenotypes"]])
}

getBioMatrixSampleWidth <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_Sample_Width"]])
}

getBioMatrixSampleHeight <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_Sample_Height"]])
}

getBioMatrixColumnPadding <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_Column_Padding"]])
}

getBioMatrixRowPadding <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_Row_Padding"]])
}

getBioMatrixGeneLabelWidth <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_Label_Width"]])
}

getBioMatrixRemarkWidth <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_Remark_Width"]])
}

getBioMatrixSummaryWidth <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_Summary_Width"]])
}

getBioMatrixDataAreaWidth<-function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_DataArea_Width"]])
}

getBioMatrixBasePositions<-function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_Base_Positions"]])
}

getBioMatrixSampleIDHeight<-function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_SampleID_Height"]])
}

getBioMatrixLegendHeight<-function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioMatrix_Legend_Height"]])
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Calculate total height of bioMatrixPlot area. 
#
#    Argument:          None
#    Returned value:    Total height of plot area in inch.
#
#    Example:   areaWidth <- getBioMatrixPlotAreaHeigth()
#
#    Last revised on September 6, 2013
#

getBioMatrixPlotAreaHeigth<-function() {

    totalGenes      <- getBioMatrixGeneNumber()
    totalPhenotypes <- getBioMatrixPhenotypeNumber()

    sampleHeight <- getBioMatrixSampleHeight()
    rowPadding   <- getBioMatrixRowPadding()
    nameHeight   <- getBioMatrixSampleIDHeight()
    legendHeight <- getBioMatrixLegendHeight()

    rowHeight      <- sampleHeight + rowPadding
    plotAreaHeight <- (totalGenes+totalPhenotypes)*rowHeight + 
                            nameHeight + legendHeight

    return (plotAreaHeight)
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Calculate total width of bioMatrixPlot area. 
#
#    Argument:  None
#    Returned value:    Total width of plot area in inch.
#
#    Example:   figureWidth <- getBioMatrixPlotAreaWidth()
#
#    Last revised on November 10, 20143
#

getBioMatrixPlotAreaWidth<-function() {

    leftNameWidth <- getBioMatrixGeneLabelWidth()
    dataAreaWidth <- getBioMatrixDataAreaWidth()
    remarkWidth    <- getBioMatrixRemarkWidth()
    summaryWidth  <- getBioMatrixSummaryWidth()
    plotAreaWidth <- leftNameWidth + dataAreaWidth + 
                        remarkWidth + summaryWidth

    return (plotAreaWidth)
}






#    __________________________________________________________________________
#    **************************************************************************
#
#   Calculate the top y coordinate for a row.
#
#   Arguments:
#       rowNumber:  non-negative integer, which row to plot 
#       area:       character vector, which area to plot
#
#   Returned value: integer, top location of a row to plot
#
#   Example: rowTop <- getBioMatrixDataRowTop(1, areaName="omicsData")
#
#   Last revised on September 5, 2014
#

getBioMatrixDataRowTop <- function(rowNumber, 
                areaName=c("omicsData", "phenotype")) {

    if(is.numeric(rowNumber) == FALSE || rowNumber<0)
        stop("Row number must be non-negative numeric.")

    areaName <- tolower(areaName)
    if(!areaName %in% c("omicsdata","phenotype"))
        stop("Incorrect plot area name.") 

    plotAreaTop     <- getBioMatrixPlotAreaHeigth()
    nameHeight      <- getBioMatrixSampleIDHeight()
    totalPhenotypes <- getBioMatrixPhenotypeNumber()

    sampleHeight <- getBioMatrixSampleHeight()
    rowPadding   <- getBioMatrixRowPadding()
    rowHeight    <- sampleHeight+rowPadding

    samplePositions <- getBioMatrixBasePositions()

    if(areaName == "phenotype") {   
        skipHeight <- (rowNumber-1)*rowHeight
    } else {  
        totalRows  <- totalPhenotypes+rowNumber-1
        skipHeight <- totalRows*rowHeight + rowPadding*2 
    }

    yStart  <- plotAreaTop-nameHeight-rowPadding-skipHeight

    return (yStart)
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Plot sample names on the top of biomatrixplot area
#
#    Argument:
#
#        sampleNames:  character vector, sample names to be plotted
#        sampleColors: character vector or R color name(s) for text color(s)
#
#    Returned value:    None
#
#    Example:   sampleNames <- colnames(exprData)
#             plotBioMatrixSampleNames(sampleNames)
#
#    Last revised on September 5, 2014
#

plotBioMatrixSampleNames <- function(sampleNames, sampleColors)
{
    totalSamples <- getBioMatrixSampleNumber()
    nameHeight   <- getBioMatrixSampleIDHeight()
    sampleWidth  <- getBioMatrixSampleWidth()
    plotHeight   <- getBioMatrixPlotAreaHeigth()

    samplePositions <- getBioMatrixBasePositions()

    xLeft <- samplePositions[,1]
    yTop  <- rep(plotHeight-nameHeight, totalSamples)

    text(xLeft-sampleWidth/2, yTop, sampleNames, pos=4, 
                srt=90, col=sampleColors)
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Plot row names on the left or right side of biomatrix plot area. Penotype
#    names, gene names, and remarks are all plotted with this function.
#
#    Argument: 
#
#        rowNames:  character vector, row names to be plotted
#        areaName:  character vector, either "omicsData" or "phenotype"
#        colors:      vector of character or R color  names
#        side:        character vector, either "left" or "right"
#        skipPlotRow: non-negative integer, total rows on plot area that 
#                     should be skipped, default 0 (means start plot from the
#                     first row)
#        skipColumns: non-negative integer, columns (sampleWidth) will be 
#                     skipped when items on remark area
#
#    Returned value:    None
#
#    Example:   geneNames <- rownames(exprData)
#               textColors <- rep("black", length(geneNames))
#               plotBioMatrixRowNames(sampleNames, colors=textColors)
#
#   Last revised on September 5, 2014
#

plotBioMatrixRowNames <- function(geneNames, areaName, colors, side="left",
                skipPlotRows=0, skipPlotColumns=0) {

    areaName <- tolower(areaName)
    if(!areaName %in% c("omicsdata","phenotype"))
        stop("Incorrect plot area name.") 

    side <- tolower(side)
    if(!side %in% c("left","right")) stop("Incorrect side definition.")

    totalGenes   <- getBioMatrixGeneNumber()
    sampleHeight <- getBioMatrixSampleHeight()

    rowNameX <- getBioMatrixGeneLabelWidth()
    position <- 2 

    if(side == "right") {
        sampleWidth   <- getBioMatrixSampleWidth()
        dataAreaWidth <- getBioMatrixDataAreaWidth()
        rowNameX <- rowNameX + dataAreaWidth + sampleWidth*skipPlotColumns
        position <- 4
    }

    for(aRow in seq_along(geneNames)) {

        yStart <- getBioMatrixDataRowTop(aRow+skipPlotRows, 
                        areaName=areaName)
        text(rowNameX,  yStart-sampleHeight/2, geneNames[aRow], 
                        pos=position, col=colors[aRow])
    }
}






#    __________________________________________________________________________
#    **************************************************************************
#
#    Plot one row of rectangles. This could be used for phenotype plot , 
#    (different background), sample layout (same background), ans heatmap 
#    (different background). 
#
#    Arguments:
#
#        rowNumber:    non-negative integer, index of the row to be plotted
#        areaName:     character vector, either "phenotype" or "omicsdata"
#        fillColor:    character vector for color names or vector of R color 
#                      specification for rectangle fill
#        borderColor:  character vector or R colors specification for border
#        topAdjust:   non-negative numeric, height that will be reduced 
#                      from row top
#        bottomAdjust: non-negative numeric, height that will be reduced 
#                      from row bottom
#
#    Return value:  None
#    Example:   plotBioMatrixSampleData(rowNumber=1, areaName="phenotype")
#
#    Last revised on September 5, 2014
#
#

plotBioMatrixSampleData <- function(rowNumber, areaName, fillColor=NA, 
                borderColor=NA, topAdjust=0, bottomAdjust=0) {

    if(is.numeric(rowNumber) == FALSE || is.numeric(topAdjust) == FALSE || 
        is.numeric(bottomAdjust) == FALSE || rowNumber<0 || topAdjust<0 || 
        bottomAdjust<0) {  
        stop("RowNumber/topAdjust/bottomAdjust must be non-negative numeric.")
    }

    areaName <- tolower(areaName)
    if(!areaName %in% c("omicsdata","phenotype"))
        stop("Incorrect plot area name.") 

    samplePositions <- getBioMatrixBasePositions()
    sampleHeight    <- getBioMatrixSampleHeight()
    totalSamples    <- getBioMatrixSampleNumber()

    yStart  <- getBioMatrixDataRowTop(rowNumber, areaName)

    xLeft   <- samplePositions[, 1]
    xRight  <- samplePositions[, 3]
    yTop    <- rep(yStart, totalSamples) - topAdjust
    yBottom <- yTop - sampleHeight + topAdjust + bottomAdjust

    rect(xLeft, yBottom, xRight, yTop, col=fillColor, border=borderColor)
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Show biomatrix layout before plot data values.
#
#    Arguments:
#
#    geneNames:     character vector, row names at most left of each data row 
#    sampleNames:   character vector, column names at top of each data 
#                   column
#    phenotypes:    character vector, name(s) of each phenotype row
#    sampleColors:  character vector of color names or vector of R color 
#    geneColors:    character vector of color names or vector of R color 
#    phenoColors:   character vector of color names or vector of R color 
#
#    Returned value:    None
#
#    Example:    geneNames <- rownames(exprData)
#                sampleNames <- colnames(exprData)
#                showBioMatrixPlotLayout(geneNames, sampleNames, "Tissue type")
#
#    Last updated on September 4, 2014
#

showBioMatrixPlotLayout <- function(geneNames, sampleNames, phenotypes,
        secondGeneNames=NULL, sampleColors=NULL, 
        geneColors=NULL, phenoColors=NULL) {  

    if(is.character(geneNames) == FALSE || 
        is.character(sampleNames) == FALSE ||
        is.character(phenotypes) == FALSE)
            stop("Arguments must be character vector.")
    
    totalGenes      <- getBioMatrixGeneNumber()
    totalSamples    <- getBioMatrixSampleNumber()
    totalPhenotypes <- getBioMatrixPhenotypeNumber()

    if(length(geneNames)!= totalGenes || 
            length(sampleNames)!=totalSamples ||
            length(phenotypes)!= totalPhenotypes)
        stop("Incorrect number of gene names/sample names/phenotypes.")

    if(is.null(secondGeneNames)==FALSE) { 
        if(is.character(secondGeneNames) == FALSE )
            stop("secondGeneNames must be character vector.")
        if(length(secondGeneNames) != totalGenes)
            stop("Length of gene names and second gene names differ.")
    }

    setBioMatrixPlotArea()

    if(is.null(sampleColors)) sampleColors <- rep("black", totalSamples)
    plotBioMatrixSampleNames(sampleNames, sampleColors)

    if(is.null(phenoColors)) phenoColors <- rep("black", totalPhenotypes)
    plotBioMatrixRowNames(phenotypes, "phenotype", phenoColors, "left")

    if(is.null(geneColors)) geneColors <- rep("black", totalGenes)
    plotBioMatrixRowNames(geneNames, "omicsData", geneColors, "left")

    if(is.null(secondGeneNames)==FALSE)
           plotBioMatrixRowNames(secondGeneNames, "omicsData", geneColors, 
                side="right", skipPlotColumns=0)

    sampleColors <- rep("lightskyblue", totalSamples)
    for(aPheno in seq_along(phenotypes)) {
        plotBioMatrixSampleData(aPheno, areaName="phenotype", 
            fillColor=sampleColors) 
    }

    sampleColors <- rep("lightcyan", totalSamples)
    for(aGene in seq_along(geneNames)) {

        plotBioMatrixSampleData(aGene, areaName="omicsData", 
                fillColor=sampleColors)
    }
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Plot heatmaps with expression data. Color scales must be user defined or  
#    based on whole dataset (NULL for max and min values). This plot all rows 
#    of input data. For single row, make color vector and call function
#    plotBioMatrixSampleData()
#
#    Arguments:
#
#    exprData:     numeric matrix (log2 values)
#    topAdjust:    non-negative numeric, height of top y coordinate should 
#                  be reduced to show different layers, default 0
#    bottomAdjust: non-negative numeric, height of bottom y coordinate should 
#                  be reduced for a small rectangle, default 0
#    maxValue:     numeric, maximum value for highest color in heatmap, set  
#                  to NULL to use the maximum value in expression dataset
#    minValue:     numeric, minimum value for lowest color in heatmap, set to  
#                  NULL to use the minimum value in expression dataset
#    heatmapColor: character vector, either "BlueWhiteRed", "GreenWhiteRed", 
#                  "GreenYellowRed", "GreenBlackRed" , or "YellowToRed"
#    skipPlotRow:  non-negative integer, total rows on plot area that should 
#                  be skipped, default 0 (start plot from the first row)
#
#    Returned value:    None
#
#    Example:   plotBioMatrixHeatmap(log2(exprValues))
#
#    Last revised on September 10, 2014
#

plotBioMatrixHeatmap <- function(exprData, topAdjust=0, bottomAdjust=0, 
            maxValue=NULL, minValue=NULL, heatmapColor="BlueWhiteRed", 
            skipPlotRow=0) {

    errMSG <- "Heatmap colors should be one from: ,
 
        BlueWhiteRed,  
        GreenWhiteRed,  
        GreenYellowRed, 
        GreenBlackRed,  
        YellowToRed. 

        Please redefine the plotColors."

    if(length(heatmapColor)>1)  stop(errMSG)

    totalGenes <- getBioMatrixGeneNumber()
    if(skipPlotRow<0 || skipPlotRow>=totalGenes) 
            stop("Incorrect row number to be skipped.")

    if(is.null(minValue) == FALSE && is.null(maxValue) == FALSE) { 
        dataRange <- c(minValue, maxValue)  
    } else {
        dataRange <- c(min(exprData), max(exprData))
    }

    colorRamp <- getHeatmapColorScales(heatmapColor)
    colorLevel <- seq(dataRange[1], dataRange[2], length=length(colorRamp))

    sampleColors <- rep("white", ncol(exprData))
    for(aRow in seq_len(nrow(exprData))) {

         for(aSample in seq_along(sampleColors)) {

            if(is.na(exprData[aRow, aSample])) {
                sampleColors[aSample] <- "grey"; next 
            }

            the.level <- which(colorLevel>=exprData[aRow, aSample])
            sampleColors[aSample] <- colorRamp[min(the.level)]
        }

        plotBioMatrixSampleData(aRow+skipPlotRow, areaName="omicsData", 
            sampleColors, topAdjust=topAdjust, bottomAdjust=bottomAdjust)
    }
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    Draw rectangle outline for one row of samples to represent categorical 
#    values. This function highlight all samples on each row.
#
#    Argument:
#
#        categoryData: vector or matrix of categorical values, e.g., "High", 
#                      '"low", and "No"
#        areaName:     character vector, either "omicsData", or "phenotype"
#        sampleColors: character vector for color names or vector of R color 
#                      specification, if defined, its length must be same as
#                      number of categories
#        lineWidth:    graphic parameter for lwd (line width), default 1
#        skipPlotRow:  non-negative integer, total rows on plot area that 
#                      should be skipped, default 0 (start from the first row)
#
#    Returned value:    None
#
#    Example:   plotBioMatrixCategoryData(categoryData, areaName="omicsData")
#
#    Last revised on September 10, 2014
#

plotBioMatrixCategoryData <- function(categoryData, 
                    areaName=c("omicsData", "phenotype"), 
                    sampleColors=palette(), lineWidth=1, skipPlotRow=0) {

    if(is.numeric(lineWidth) == FALSE || lineWidth<0)
        stop("LineWidth must be non-negatice nemeric.") 

    areaName <- tolower(areaName)
    if(!areaName %in% c("omicsdata","phenotype"))
        stop("Incorrect plot area name.") 

    categoryNames <- unique(as.vector(categoryData))
    categoryNames <- sort(categoryNames)

    if(length(which(is.na(categoryNames)))>0)
        categoryNames   <- categoryNames[-which(is.na(categoryNames))]
    totalcategories <- length(categoryNames)

    if(length(sampleColors)<totalcategories)
        stop("Length of sample colors is less than total categories.")

    areaName <- tolower(areaName)
    if(!areaName %in% c("omicsdata","phenotype"))
        stop("Incorrect plot area name.") 

    samplePositions <- getBioMatrixBasePositions()
    sampleHeight    <- getBioMatrixSampleHeight()

    for(aRow in seq_len(nrow(categoryData))) {
        yStart  <- getBioMatrixDataRowTop(aRow+skipPlotRow, areaName)

        xLeft   <- samplePositions[, 1]
        xRight  <- samplePositions[, 3]
        yTop    <- rep(yStart, length(xLeft))
        yBottom <- yTop - sampleHeight

        plotData <- categoryData[aRow, ]
        for(aCategory in seq_len(totalcategories) ) {
            sampleIndex <- which(plotData == categoryNames[aCategory])
            rect(xLeft[sampleIndex], yBottom[sampleIndex], 
                xRight[sampleIndex], yTop[sampleIndex], 
                col=NA, border=sampleColors[aCategory], lwd=lineWidth)

            if(lineWidth>1)
                rect(xLeft[sampleIndex], yBottom[sampleIndex], 
                    xRight[sampleIndex], yTop[sampleIndex], 
                    col=NA, border="white", lwd=1)
        }
    }
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    Plot binary data as points in the inside of each rectangle(sample). The
#    point type and colors can be defined by user. This function plot all rows
#    on omics data area and only the positive samples will be shown with 
#    colored points. For one row plot, pass data as vector and supply correct
#    skipPlotRow parameter to define where to plot.
#
#    Arguments:
#
#        binaryData:  integer vector or matrix with 0/1 only 
#        areaName:    character vector, either "omicsData" or "phenotype" 
#        scatterType: non-negative integer, same as pch, default 19 
#        scatterSize: non-negative numeric, same as cex 
#        totalSubRow: non-negative integer, how many sub-rows in a sample area
#        subRowIndex: non-negative integer, which subrow will be plotted 
#        sampleColor: character vector for color name or R color specification
#        skipPlotRow: non-negative integer, total rows on plot area that 
#                     should be skipped, default 0 (start from the first row)
#
#    Returned value:    None
#
#    Example:   plotBioMatrixBinaryData(binaryData, areaName="omicsData")
#
#    Last revised on September 10, 2014
#
#

plotBioMatrixBinaryData <- function(binaryData, areaName="omicsData", 
        scatterType=19, scatterSize=1, totalSubRow=1, subRowIndex=1, 
        sampleColor="black", skipPlotRow=0) {

    binaryValues <- unique(as.vector(binaryData))
    if(length(binaryValues)!=2) stop("Input data is not binary.")

    if(!scatterType %in% 1:25 ) stop("Incorrect scatter type.") 
    if(!scatterSize %in% 1:3)   stop("Incorrect scatter size.") 

    if(totalSubRow>4 || totalSubRow<1) stop("Incorrect sub row size.")
    if(subRowIndex<1 || subRowIndex>totalSubRow) 
        stop("Incorrect sub row location.") 

    totalGenes <- getBioMatrixGeneNumber()
    if(skipPlotRow>=totalGenes) 
        stop("Incorrect row number to be skipped.")

    samplePositions <- getBioMatrixBasePositions()
    sampleHeight    <- getBioMatrixSampleHeight()

    subRowHeight <- sampleHeight/totalSubRow
    pointHeight <- subRowHeight*subRowIndex - subRowHeight/2

    pointX <- (samplePositions[, 1] + samplePositions[, 3]) / 2
    for(aRow in seq_len(nrow(binaryData)))  {
        rowTop  <- getBioMatrixDataRowTop(aRow+skipPlotRow, areaName)
        pointY <- rep(rowTop - pointHeight, length(pointX))

        sampleIndex <- which(binaryData[aRow,] == 1)
        points(pointX[sampleIndex], pointY[sampleIndex], pch=scatterType, 
                    cex=scatterSize, col=sampleColor[1])
    }
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    Plot bars with frequency data such as RNASeq read coverage. This function
#    plots all rows on omics data area or one column on remark area. For single
#    row plot, pass data as vector and supply desired skipPlotRow parameter to
#    define where to plot
#
#    Arguments:
#
#        barData:   numeric matrix or vector with values in range of 0 ~ 1
#        barColor:  character vector for color name or R color specification
#        areaName:  character vector, name of plot area, currently use
#                   "omicsData" only
#        byRow:     logic, whether plot bars for each row or not
#        skipPlotRow: non-negative integer, how many row(s) to be skipped from
#                       the first row
#        skipColumns: non-negative integer, how many column(s) to be skipped
#                       from the first column
#
#    Returned value:    None
#
#    Example:   plotBioMatrixBars(barData, areaName="omicsData")
#
#    Last revised on September 10, 2014
#

plotBioMatrixBars <- function(barData, barColor="red", areaName="omicsData",
                byRow=TRUE, skipPlotRow=0, skipPlotColumns=0) {

    if(is.numeric(barData) == FALSE || max(barData)>1 || min(barData)<0)
            stop("Input data must be numeric in range of 0 ~ 1.") 

    if(byRow == FALSE) { 
        if(is.vector(barData) == FALSE) stop("Input data should be a vector.")

        totalRows  <- length(barData) 
        barsPerRow <- 1

        if(totalRows != getBioMatrixGeneNumber()) 
            stop("Incorrect input data length.")
    } else {
        if(is.vector(barData)) {
            totalRows <- 1 
            barsPerRow <- length(barData)
            barData <- matrix(barData, ncol=barsPerRow)
        } else if(is.matrix(barData)) { 
            totalRows <- nrow(barData)
            barsPerRow <- ncol(barData)
        } else { stop("Bar data should be numeric vector or matrix") }
    }

    samplePositions <- getBioMatrixBasePositions()
    sampleHeight    <- getBioMatrixSampleHeight()
    barData <- barData*sampleHeight

    barXLeft  <- samplePositions[,1]
    barXRight <- samplePositions[,3]

    if(byRow == FALSE) {
        rowNameX <- getBioMatrixGeneLabelWidth()
        dataAreaWidth <- getBioMatrixDataAreaWidth()
        sampleWidth   <- getBioMatrixSampleWidth()
        samplePadding <- getBioMatrixColumnPadding()
        columnWidth   <- sampleWidth+samplePadding

        barXLeft <- rowNameX + dataAreaWidth
        barXLeft <- barXLeft + columnWidth*skipPlotColumns
        barXRight <- barXLeft + sampleWidth
    }


    for(aRow in seq_len(totalRows)) {

        yStart  <- getBioMatrixDataRowTop(aRow+skipPlotRow, areaName)
        yBottom <- yStart - sampleHeight

        if(byRow == FALSE) { #  only one bar for this row
            rect(barXLeft, yBottom, barXRight, 
                yBottom+barData[aRow], col= barColor)
        } else {
            for(aBar in seq_len(barsPerRow)) {
                rect(barXLeft[aBar], yBottom, barXRight[aBar],
                    yBottom+barData[aRow, aBar], col=barColor)
            }
        }
    }
}

#   Last revised on May 18, 2015
#   ______________________________________________________________________
