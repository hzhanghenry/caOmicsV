#
#    R Package for Cancer Genomic Data Visualization (caOmicsV) 
#
# 
#    __________________________________________________________________________
#
#    Data types to  be visualized included:
#
#    1). Clinical information
#    2). Gene expression from microarray or Next Generation Sequencing
#    3). Copy number variations (deletions, insertions, and amplifications)
#    4). Mutations (mutation status presented with binary 0/1)
#    5). Methylations
#
#    Layouts included:
#
#    1)    BioNetCircos:    circos plot on network layout 
#    2)    Biomatix:        matrix layout
#
#    Created on August 8, 2014
#    Revised on March 18, 2015 in compliance with Bioconductor coding style
# 
#    by Hongen Zhang, Ph.D. (hzhang@mail.nih.gov)
#
#    Genetics Branch
#    Center for Cancer Research 
#    National Cancer Institute
#    National Institutes of Health
#    Bethesda, Maryland 20892
#




#    Private Environment to hold caOmicsV objects such as node diameters, 
#    coordinates of each node center, coordinates of each data point ...
#    __________________________________________________________________________
#    **************************************************************************

CA_OMICS_ENV  <- new.env()
CA_OMICS_NAME <- "CA_OMICS_ENV"
CA_OMICS_NA_STRING <- NULL





#    __________________________________________________________________________
#    **************************************************************************
#
#    A easy way to plot cancer genomic data with bioNetCircos layout
#
#    Arguments:
#
#    dataSet:        A list of objects containing all plot datasets, include:
#
#    sampleNames:    character vector, sample names to be plotted
#    geneNames:      character vector, names of genes for each row
#
#    sampleInfo:     data frame, sample information such as tissue type, 
#                    diagnosis, gender, age, et al.
#
#    heatmapData:    list of data frames, numeric data, maximum number
#                    of datasets included: 2, values should be log2
#                    transformed. The first column in each dataset must 
#                    be gene names same as geneNames above and column 
#                    names must be same as sampleNames above
#
#    categoryData:   list of data frames with categorical data, maximum
#                    number of datasets included: 2. The first column of 
#                    each dataset must be gene names same as geneNames 
#                    above and column names must be same as sampleNames
#
#    binaryData:     list of data frames with binary data (0/1), maximum
#                    number of datasets included: 3. The first column in  
#                    each dataset must be gene names same as geneNames  
#                    ablove and column names must be same as sampleNames
#
#    summaryData:    list of data frames, summarization for genes or for
#                    samples. Summary data for gene must have length of 
#                    sampleNames and summary data for samples must have 
#                    length of geneNames
#
#    graph:          an igraph object for customized graph plot. Total 
#                    number of nodes cannot be less than total number of 
#                    genes in each dataset.
#
#    heatmapColor:   character vector, one of
#
#                    BlueWhiteRed:   from blue to white then red
#                    GreenWhiteRed:  from green to white then red
#                    GreenYellowRed: from green to yellow then red
#                    GreenBlackRed:  from green to black then red
#                    YellowToRed:     from yellow to red
#                    BlackOnly:      default black colors
#
#    Returned value:    None
#
#    example:             plotBioNetCircos(dataSet)
#
#    Last revised on January 26, 2015
#

plotBioNetCircos<-function(dataSet, graph=NULL, heatmapMax=NULL, 
            heatmapMin=NULL, heatmapColor="BlueWhiteRed") {

    sampleNames  <- dataSet$sampleNames
    geneNames    <- dataSet$geneNames
    numOfSamples <- length(sampleNames)

    widthOfSample    <- 100
    widthBetweenNode <- 3
    lengthOfRadius   <- 10;

    numOfSampleInfo <- nrow(dataSet$sampleInfo) - 1
    numOfSummary <- ifelse(dataSet$summaryByRow, 0, col(dataSet$summaryInfo)-1)
    numOfHeatmap <- length(dataSet$heatmapData)
    numOfCategory <- length(dataSet$categoryData)
    numOfBinary <- length(dataSet$binaryData)

    dataNum <- sum(numOfSampleInfo, numOfSummary, numOfHeatmap, 
        numOfCategory, numOfBinary)

    trackheight <- 1.5
    widthOfPlotArea  <- dataNum*2*trackheight

    if(is.igraph(graph))  { 
            bioNet <- graph;
    } else { 
        if(numOfHeatmap<1) {stop("Data for graph not found.") }
        expr <- dataSet$heatmapData[[1]]
        bioNet <- bc3net(expr) 
    }

    initializeBioNetCircos(bioNet, numOfSamples, widthOfSample, 
            lengthOfRadius, widthBetweenNode, widthOfPlotArea)
    caOmicsVColors <- getCaOmicsVColors()
    supportedType  <- getCaOmicsVPlotTypes()

    showBioNetNodesLayout();
    par(cex=0.6);
    labelBioNetNodeNames(nodeList=1:length(geneNames),
        labelColor=caOmicsVColors[3], labelLocation="bottom", 
        labelOffset=0.75)

    eraseBioNetNode()
    inner <- lengthOfRadius/2
    outer <- inner + trackheight

    # sample data is in a data frame. First row is sample ID
    #
    for(aSam in seq_len(nrow(dataSet$sampleInfo))[-1]) {

        print(paste("plot", supportedType[1]))
        groupInfo  <- as.character(dataSet$sampleInfo[aSam, ])
        sampleType <- unique(groupInfo)

        colorSet <- caOmicsVColors
        if(length(sampleType)>length(colorSet))
                colorSet <- rainbow(length(sampleType))
        sampleColors <- rep(colorSet[1], length(sampleNames))

        for(colorItem in seq_len(length(sampleType))[-1]) {
            smapleIndex <- which(groupInfo == sampleType[colorItem])
            sampleColors[smapleIndex] <- colorSet[colorItem]
        }

        groupInfo <- matrix(groupInfo, nrow=1)
        bioNetCircosPlot(dataValues=groupInfo, supportedType[1], 
                    outer, inner, sampleColors)

        inner <- outer + 0.5;  
        outer <- inner +  trackheight  
    }

    for(aMap in seq_along(dataSet$heatmapData)) {
        print(paste("plot", supportedType[4]))

        exprData <- dataSet$heatmapData[[aMap]]
        bioNetCircosPlot(exprData, supportedType[4], outer, inner, 
                plotColors=heatmapColor, heatmapMax, heatmapMin)
        inner <- outer + 0.5;
        outer <- inner +  trackheight  
    }

    for(aGroup in seq_along(dataSet$categoryData)) {

        print(paste("plot", supportedType[2]))

        categoryData <- dataSet$categoryData[[aGroup]]
        plotColors <- rep(caOmicsVColors[1], ncol(categoryData))
        bioNetCircosPlot(categoryData, supportedType[2], outer, 
            inner, plotColors)

        inner <- outer + 0.5;
        outer <- inner +  trackheight;  
    }

    for(aGroup in seq_along(dataSet$binaryData)) {

        print(paste("plot", supportedType[3]))

        binaryData <- dataSet$binaryData[[aGroup]]
        plotColors <- rep(caOmicsVColors[1], ncol(binaryData))
        bioNetCircosPlot(binaryData, supportedType[3], outer, 
            inner, plotColors)

        inner <- outer + 0.5;
        outer <- inner +  trackheight  
    }
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    A easy way to plot cancer genomic data with bioMAtrix layout
#
#    Arguments:
#
#    dataSet:       A list object containing all plot datasets, include:
#
#    sampleNames:    character vector, sample names to be plotted
#    geneNames:      character vector, names of genes for each row
#
#    sampleInfo:     data frame, sample information such as tissue type, 
#                    diagnosis, gender, age, et al.
#
#    heatmapData:    list of data frames, numeric data, maximum number
#                    of datasets included: 2, values should be log2
#                    transformed. The first column in each dataset must 
#                    be gene names same as geneNames above and column 
#                    names must be same as sampleNames above
#
#    categoryData:   list of data frames with categorical data, maximum
#                    number of datasets included: 2. The first column of 
#                    each dataset must be gene names same as geneNames 
#                    above and column names must be same as sampleNames
#
#    binaryData:     list of data frames with binary data (0/1), maximum
#                    number of datasets included: 3. The first column in  
#                    each dataset must be gene names same as geneNames  
#                    ablove and column names must be same as sampleNames
#
#    summaryData:    list of data frames, summarization for genes or for
#                    samples. Summary data for gene must have length of 
#                    sampleNames and summary data for samples must have 
#                    length of geneNames
#
#    graph:          an igraph object for customized graph plot. Total 
#                    number of nodes cannot be less than total number of 
#                    genes in each dataset.
#
#    heatmapColor:   character vector, one of
#
#                   BlueWhiteRed:   from blue to white then red
#                   GreenWhiteRed:  from green to white then red
#                   GreenYellowRed: from green to yellow then red
#                   GreenBlackRed:  from green to black then red
#                   YellowToRed:    from yellow to red
#                   BlackOnly:      default black colors
#
#
#    Returned value: None
#    Example:   plotBioMatrix(dataSet, summaryType="text", summarybyRow=TRUE);
#
#   Last revised on November 24, 2014
#

plotBioMatrix <- function(dataSet, summaryType=c("text", "bar"), 
            summarybyRow=TRUE, heatmapMax=NULL, heatmapMin=NULL, 
            heatmapColor="BlueWhiteRed") {

    if(is.null(dataSet$geneNames)) stop("Gene names must be provided.")
    numOfGenes <- length(dataSet$geneNames);

    if(is.null(dataSet$sampleNames)) stop("Sample names must be provided.")
    numOfSamples <- length(dataSet$sampleNames); 

    numOfPhenotypes <- nrow(dataSet$sampleInfo)-1;
    if(numOfPhenotypes<1) stop("Sample info must have two or more columns.")

    numOfHeatmap    <- length(dataSet$heatmapData);
    if(numOfHeatmap>2) stop("Number of headmap data is limited to 2.")

    numOfSummary    <- length(dataSet$summaryData);
    phenotypes      <- rownames(dataSet$sampleInfo)[-1];

    sampleHeight=0.4;
    sampleWidth=0.1; 
    samplePadding=0.025;
    geneNameWidth=1;
    sampleNameHeight=5;
    remarkWidth=2; 
    summaryWidth=1;
    rowPadding=0.1;

    initializeBioMatrixPlot(numOfGenes, numOfSamples, numOfPhenotypes, 
        sampleHeight, sampleWidth, samplePadding,  rowPadding, 
        geneNameWidth, remarkWidth, summaryWidth, sampleNameHeight);

    par(cex=0.75);
    showBioMatrixPlotLayout(dataSet$geneNames, dataSet$sampleNames, 
        phenotypes, dataSet$secondGeneNames);
    caOmicsVColors <- getCaOmicsVColors()

    for(aType in seq_len(numOfPhenotypes)) {

        rowIndex <- aType+1;
        sampleGroup <- as.character(dataSet$sampleInfo[rowIndex,])
        sampleTypes <- unique(sampleGroup);

        plotColors <- caOmicsVColors
        if (length(sampleTypes)>length(plotColors)) 
            plotColors <- rainbow(length(sampleTypes))

        sampleColors <- rep(plotColors[1], length(sampleGroup))
        for(aColor in seq_len(length(sampleTypes))[-1]) {  
            sampleIndex <- grep(sampleTypes[aColor], sampleGroup)
            sampleColors[sampleIndex] <- plotColors[aColor]   
        }

        plotBioMatrixSampleData(aType, "phenotype", sampleColors)

        geneLabelX <- getBioMatrixGeneLabelWidth()
        maxAreaX   <- getBioMatrixDataAreaWidth()
        legendH    <- getBioMatrixLegendHeight()
        plotAreaH  <- getBioMatrixPlotAreaHeigth()
        sampleH    <- getBioMatrixSampleHeight()

        sampleLegendX <- geneLabelX + maxAreaX 
        sampleLegendY <- plotAreaH + legendH - length(sampleTypes)*sampleH
        colors <- plotColors[1:length(sampleTypes)]
        legend(sampleLegendX, sampleLegendY, legend=sampleTypes, 
                    fill=colors,  bty="n", xjust=0)
    }

    for(aHeatmap in seq_len(numOfHeatmap)) {

        topStart <- ifelse(aHeatmap==1, 0, sampleHeight/2)

        heatmapData <- dataSet$heatmapData[[aHeatmap]]
        plotBioMatrixHeatmap(heatmapData, topAdjust=topStart,
                    maxValue=heatmapMax, minValue=heatmapMin)
    }

    for(aCat in seq_along(dataSet$categoryData)) {
        categoryData <- dataSet$categoryData[[aCat]]
        totalCategory <- length(unique(as.numeric(categoryData)))

        plotColors <- rev(getCaOmicsVColors())
        if(totalCategory > length(plotColors)) { 
            stop("Too many categories to plot.") 
        }

        plotBioMatrixCategoryData(categoryData, areaName="omicsData", 
            sampleColors=plotColors[1:totalCategory])
    }

    for(aBinary in seq_along(dataSet$binaryData)) {

        binaryData <- dataSet$binaryData[[aBinary]];
        plotBioMatrixBinaryData(binaryData, sampleColor=caOmicsVColors[4]);
    
        if(length(dataSet$binaryData)>1) {
            binaryData <- dataSet$binaryData[[2]];
            binaryData <- as.matrix(binaryData[,-1])
            plotBioMatrixBinaryData(binaryData, sampleColor=caOmicsVColors[3])
        }
    }

    for(aSum in seq_along(dataSet$summaryInfo)) {

        summaryData  <- dataSet$summaryInfo[[aSum]][, 2];
        summaryTitle <- colnames(dataSet$summaryInfo[[aSum]])[2];

        if(length(dataSet$heatmapData)>1) {
                remarkWidth <- getBioMatrixRemarkWidth();
                sampleWidth <- getBioMatrixSampleWidth();
                col2skip <- remarkWidth/sampleWidth + 2;
        } else {  col2skip <- 1; }

        plotBioMatrixRowNames(summaryTitle, areaName="phenotype", 
            colors="black", side="right", skipPlotColumns=col2skip);
                
        if(summaryType == "text") {
            plotBioMatrixRowNames(summaryData, "omicsData", 
                colors=caOmicsVColors[3], side="right", 
                skipPlotColumns=col2skip)
        }  else {
            plotBioMatrixBars(summaryData, caOmicsVColors[1], 
                areaName="omicsData", byRow=TRUE);
        }
    }
}



#   ________________________________________________________________________
#   ************************************************************************
#
#   Methods to get and set NA strings used by caOmicsV package
#

getDefaultNaStrings <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["CA_OMICS_NA_STRING"]])
}

setDefaultNaStrings <- function(nullStrings) {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())
    caOmicsVEnvironment[["CA_OMICS_NA_STRING"]] <- nullStrings
}






#   ________________________________________________________________________
#   ************************************************************************
#
#   Methods to get and set default plot colors for caOmicsV plot
#

getCaOmicsVColors <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return(caOmicsVEnvironment[["CA_OMICS_COLORS"]])
}

setCaOmicsVColors <- function(colorList=NULL)
{
    if(is.null(colorList)) {  
        colorList <- c("red", "blue", "black", "green", "cyan", 
                "brown", "magenta", "gold") 
    }

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())
    caOmicsVEnvironment[["CA_OMICS_COLORS"]] <- colorList
}



#    __________________________________________________________________________
#    **************************************************************************
#
#    Method to get supported plot types
#

getCaOmicsVPlotTypes <- function() {

    supportedPlotTypes <- c("polygon", "bar", "points", "heatmap", 
        "line", "category", "binary")

    return (supportedPlotTypes)
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Make plot dataset from data frames so that the plot could be handled 
#    automatically.
#
#    The input data should included:
#
#    1.  Sample information: required
#    2.  Continue (expression) numeric data for heatmap plot. Maximum two 
#        for bioMatrix layout
#    3.  One or more categorical dataset for colored outline of samples 
#        for bioMatrix layout or colored polygons for bioNetCircos layout
#    4.  One or more (total no more than 3) of binary dataset for point 
#        plots in bioMatrix layout or colored polygons in bioNetCircos layout
#    5.  Summary data (usually percentage in number or text) for extra rows  
#        below omics data area or extra columns on remark area in biomatrix 
#        layout or bar plot in bioNetCircos layout
#
#    All dataset should be sorted by row names and column names where apply. 
#    This could be done with supplemented methods in this script.
#
#    Arguments:
#
#        sampleNames:   character vector, names of samples to plot
#        geneNames:     character vector, names of genes to plot
#        sampleData:    data frame, sample information, with rows for samples 
#                       and columns for features. The first column must be 
#                       sample names
#        heatmapData: list of data frames, continue numeric data, maximum of 2
#        categoryData: list of data frames, categorical data, maximum of 2 
#        binaryData:  list of data frames, binary data (0/1), maximum of 2
#        summaryData: list of data frames, summarization for genes or samples
#        secondGeneNames: character vector, names of second set of genes to 
#                     select
#
#   Returned value: A list containing all data objects. Omics data will 
#                      be in data matrix with row and column names
#
#   Example:    dataSet <- getPlotDataSet(sampleNames, geneNames, sampleData, 
#                            heatmapData=list(A, B), categoryData=list(C, D), 
#                            binaryData=list(E, F), summaryData=list(G, H))
#
#   Last revised on March 4, 2015
#

getPlotDataSet<-function(sampleNames, geneNames, sampleData, heatmapData=list(), 
            categoryData=list(), binaryData=list(), summaryData=list(),
            secondGeneNames=NULL ) {

    numDataset <- c(length(heatmapData),length(categoryData),
                        length(binaryData))
    if(sum(numDataset) == 0) { 
        stop("There is at least one omics data to plot.") 
    }

    if(is.null(sampleNames) || is.null(geneNames)) 
    { stop("Missing sampleNames or geneNames.") }

    sampleData <- data.frame(t(sampleData))
    sampleData <- as.matrix(sampleData)

    numOfHeatmapData <- length(heatmapData);
    if(numOfHeatmapData>0) {
        for(aData in 1:numOfHeatmapData) {
            theData <- heatmapData[[aData]]
            rownames(theData) <- as.character(theData[,1])
            heatmapData[[aData]] <- as.matrix(theData[,-1])
        }
    }

    numOfCategoryData <- length(categoryData);
    if(numOfCategoryData>0) {
        for(aData in 1:numOfCategoryData) {
            theData <- categoryData[[aData]]
            rownames(theData) <- as.character(theData[,1])
            categoryData[[aData]] <- as.matrix(theData[,-1])
        }
    }

    numOfBianryData <- length(binaryData)
    if(numOfBianryData>0) {
        for(aData in 1:numOfBianryData) {
            theData <- binaryData[[aData]]
            rownames(theData) <- as.character(theData[,1])
            binaryData[[aData]] <- as.matrix(theData[,-1])
        }
    }

    numOfSummaryData <- length(summaryData)
    if(numOfSummaryData>0) {
        for(aData in 1:numOfSummaryData) {
            theData <- summaryData[[aData]]
            rownames(theData) <- as.character(theData[,1])
            summaryData[[aData]] <- as.matrix(theData)
        }
    }

    return (list(sampleNames=sampleNames, geneNames=geneNames, 
                secondGeneNames=secondGeneNames, sampleInfo=sampleData, 
                heatmapData=heatmapData, categoryData=categoryData, 
                binaryData=binaryData, summaryInfo=summaryData) )
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Extract subset of sample information.
#
#    Arguments:
#
#        sampleNames: character vector, names of samples to select
#        sampleData:  data frame, column 1 must be sample names
#
#    Returned value:    a data frame with subset of input sample data with the 
#                       order same as sample names
#
#    Example:    sampleNames <- colnames(sampleData)[1:10]
#                sampleInfor <- getPlotSampleData(sampleData, sampleNames)
#
#    Last revised on October 24, 2014
#
#

getPlotSampleData<-function(sampleData, sampleNames) {

    if(is.data.frame(sampleData) == FALSE ) { 
        stop("sampleData must be a data frame.") 
    }

    if(length(sampleNames) == 0)  stop("Sample names missing.") 

    columns <- which(as.character(sampleData[,1]) %in% sampleNames)
    if(length(columns)!=length(sampleNames))
        stop("Missing or redundant samples found in sample data.")

    sampleData   <- sampleData[columns,]

    sampleNameOrder <- order(sampleNames);
    rowOrder <- order(as.character(sampleData[,1]))

    sampleData <- sampleData[rowOrder , ]
    sampleData <- sampleData[order(sampleNameOrder) , ]

    return (sampleData)
}



#    __________________________________________________________________________
#    **************************************************************************
#
#    Extract required rows and columns from a omics data set and sort both row 
#    and columns based on the order of gene names and sample names.
#
#    Arguments:
#
#        omicsData: data frame with all samples and all genes
#        colNames:  character vector, names of columns to be extracted
#        rowNames:  character vector, names of rows to be extracted
#
#    Returned value:    data frame with subset of input data
#
#    Example:   exprData <- getPlotData(omicsData, colNames, rowNames)
#
#    Last revised on September 12, 2014
#

getPlotOmicsData<-function(omicsData, sampleNames, geneNames) {

    if(!is.data.frame(omicsData)) stop("OmicsData must be in data frame.")

    totalSamples <- length(sampleNames)
    totalGenes <- length(geneNames)

    columns <- which(colnames(omicsData) %in% sampleNames)
    if(length(columns)!=totalSamples) stop("Unable to match all samples.")
    subSet <- omicsData[, c(1, columns)]

    rows <- which(as.character(subSet[,1]) %in% geneNames)
    if(length(rows)!=totalGenes) stop("Unable to match all genes.")
    subSet <- subSet[rows, ]

    subSet <- sortOmicsDataByColumn(subSet, sampleNames)
    subSet <- sortOmicsDataByRow(subSet, geneNames)

    return (subSet)
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Extract required rows and columns from a summary data set sorted by the 
#    order of genes names or the order of sample names
#
#    Arguments:
#
#        summaryData: data frame with summary data for each gene (rows are for
#                     genes and columns are summary values) or for each sample
#                     (rows are summary values and columns are sample names)
#        sampleNames: character vector, names of samples to be extracted
#        geneNames:  character vector, names of gene to be extracted
#
#   Returned value: data frame with subset of input data 
#
#   Example: mutatPercentage <- getPlotSummaryData(summaryData, 
#                                    sampleNames, geneNames)
#
#   Last revised on September 12, 2014
#

getPlotSummaryData <- function(summaryData, sampleNames=NULL, geneNames=NULL) {

    if(is.data.frame(summaryData) == FALSE) 
    { stop("Summary data must be in data frame.") }

    if(is.null(sampleNames) && is.null(geneNames))
    { stop("Either sample names or gene names must be defined.") }

    sampleID <- colnames(summaryData)[-1]
    geneID   <- as.character(summaryData[,1])

    if(is.null(geneNames)) {
        columns <- which(sampleID %in% sampleNames)
    if(length(columns) == 0) { stop("No column found.") }

        subSet  <- summaryData[, c(1, columns+1)]
        subSet  <- sortOmicsDataByColumn(subSet, sampleNames)

    } else if(is.null(sampleNames)) {

        rows <- which(geneID %in% geneNames)
    if(length(rows) == 0) { stop("No row found.") }

        subSet  <- summaryData[rows, ]
        subSet <- sortOmicsDataByRow(subSet, geneNames)

    } else { stop("Unable to get subset from summary data.") }


    return (subSet);
}



#    __________________________________________________________________________
#    **************************************************************************
#
#    Extract a subset of plot data that are associated to other plot data such 
#    as miRNA expression or DNA copy number variants that are related to a set 
#    of differentially expressed genes 
#
#    Arguments:
#
#    omicsData:   data frame, the dataset from which subset is extracted
#    linkData:    data frame, gene names and their related items. The first
#                 column must be the item to which the second is linked to
#    geneNames:   character vector, names of genes to which the second item
#                 are linked
#
#    Returned value:    a data frame of subset data
#
#    Example: plotdata <- getRelatedPlotData(omicsData=miRNASeq, 
#                           linkData=RNA_miRNA_link, geneNames=deGenes);
#

getRelatedPlotData<-function(omicsData, linkData, geneNames) {

    if(is.data.frame(omicsData) == FALSE || is.data.frame(linkData) == FALSE) 
        stop("OmicsData and link data must be in data frame.")

    if(is.character(geneNames) == FALSE || is.vector(geneNames) == FALSE ||
        length(geneNames) == 0) {
            stop("Incorrect gene names defined."); 
    }

    totalGenes <- length(geneNames);
    totalSams  <- ncol(omicsData)-1;
    totalCol   <-  ncol(omicsData);

    subset <- data.frame(linkData[, 2], matrix(rep(0, totalGenes*totalSams), 
                            ncol=totalSams));
    colnames(subset) <- colnames(omicsData);

    for(aRow in 1:length(geneNames)) {

        aGene <- as.character(geneNames[aRow]);
        linkRow <- which(as.character(linkData[,1]) == aGene);
        anItem <- as.character(linkData[linkRow , 2])

        omicsRow <- which(as.character(omicsData[,1]) == anItem);
        if(length(omicsRow)!=1) stop("Incorrect defined items in omicsData.")

        subset[aRow, ] <- omicsData[omicsRow, ]
    }


    return (subset);
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Sort the clinical data by one column defined by the byItem. 
#
#    Arguments:
#
#        clinicalData: A data frame with rows for samples and columns for 
#                        features. Sample names must be in the first column.
#        byItem: Character vector of a group feature (column header)
#
#    Return value: ordered clinical dataset
#
#    Example:   clinicalData <- sortClinicalData(clinicalData,"diagnosis")
#
#    Last revisited on June 19, 2014
#
#

sortClinicalData <- function(clinicalData, byItem) {

    if(is.null(clinicalData) || is.null(byItem))
        stop("Missing input data or sort item.\n") 

    if(is.data.frame(clinicalData) == FALSE)
        stop("Input dataset must be data frame!\n") 

    if(is.character(byItem) == FALSE || is.vector(byItem) == FALSE)
        stop("Sort item must be a character vector.\n") 

    theCol <- which(colnames(clinicalData) == byItem)
    if(length(theCol) == 0) 
        stop(paste(byItem, "column not find!\n"))

    if(length(theCol)>1)
        stop(paste("Find more than one column for", byItem, "!\n"))

    groupData <- clinicalData[order(as.character(clinicalData[, 1])), ]
    clinicalData <- groupData[order(groupData[, theCol]), ]

    return(clinicalData)
}



#    __________________________________________________________________________
#    **************************************************************************
#
#    Sort the omics data by column header so that the columns of omics data  
#    will be in a specific order, e.g, all tumors followed by all normals.
#
#    Arguments:
#
#        omicsData: A data frame that holds genomic data such as gene 
#                    expression, SNV, RNASeq ....  The column must be 
#                    the sample names that are same as the sample names 
#                    in clinical data.
#
#        sampleNames: Character vector, sample names in a given order 
#                     (such as diagnosis). No redundant allowed
#
#    Return value: Ordered omics dataset in a data frame
#
#    Example:   omicsData <- sortOmicsDataByColumn(omicsData, sampleNames)
#
#    Last revisited on October 20, 2014
#
#

sortOmicsDataByColumn <- function(omicsData, sampleNames) {

    if( is.null(omicsData) || length(sampleNames) == 0) 
            stop("Missing input data or sort item.\n") 

    if( is.data.frame(omicsData) == FALSE)  
            stop("Input data must be in a data frame!\n") 

    if((ncol(omicsData)-1) != length(sampleNames))  
            stop("Sample sizes different!\n") 

    if(length(unique(sampleNames)) != length(sampleNames))
            stop("Find redundant column headers!\n")

    omicsSamples <- colnames(omicsData)[-1]
    if(length(unique(omicsSamples)) != length(omicsSamples))
        stop("Find redundent sample names in omics data!\n")

    sampleNameOrder <- order(sampleNames)
    omicsSampleOrder <- order(omicsSamples);

    groupSamples <- sampleNames[sampleNameOrder ]
    omicsSamples  <- omicsSamples[omicsSampleOrder];
    if(sum(groupSamples == omicsSamples) != length(sampleNames))
        stop("Sample name mismatch Found."); 

    omicsData <- omicsData[, c(1, omicsSampleOrder+1)]
    omicsData <- omicsData[, c(1, order(sampleNameOrder)+1)]

    omicsSamples <- colnames(omicsData)[-1]
    if(sum(sampleNames == omicsSamples) != length(sampleNames)) 
        stop("This message should never be reported!\n") 

    return (omicsData)
}



#    __________________________________________________________________________
#    **************************************************************************
#
#    Sort omics data by row (genes) so that the rows of omics data is in an  
#    specific order same as the geneNames
#
#    Arguments:
#
#        omicsData: A data frame that holds genomic data such as gene 
#                   expression, SNV, RNASeq ....  The column must be the 
#                   sample names that are same as the sample names in 
#                   clinical data.
#
#        geneNames: Character vector, gene names in a given order (such as by
#                    p values). Redundant gene names is allowed.
#
#    Return value:  Ordered omics dataset in a data frame
#
#    Example:   omicsData <- sortOmicsDataByRow(omicsData, geneNames)
#
#    Last revisited on October 20, 2014
#

sortOmicsDataByRow <- function(omicsData, geneNames) {

    if( is.null(omicsData) || length(geneNames) == 0) 
        stop("Missing input data or sort item.\n") 

    if( is.data.frame(omicsData) == FALSE)  
        stop("Input data must be in a data frame!\n") 

    if(nrow(omicsData) != length(geneNames))  
        stop("Row sizes different!\n") 

    omicsGenes <- as.character(omicsData[, 1])
    geneNameOrder <- order(geneNames)
    omicsGeneOrder <- order(omicsGenes)

    theGeneNames <- geneNames[geneNameOrder]
    theomicsGenes  <- omicsGenes[omicsGeneOrder];
    if(sum(theGeneNames == theomicsGenes) != length(geneNames))
        stop("Gene name mismatch Found.")

    omicsData <- omicsData[omicsGeneOrder, ]
    omicsData <- omicsData[order(geneNameOrder), ]

    omicsGenes <- as.character(omicsData[, 1])
    if(sum(geneNames == omicsGenes) != length(geneNames)) 
        stop("This message should never be reported!\n")

    return (omicsData)
}






#    __________________________________________________________________________
#    **************************************************************************
#
#    Print out supported circos plot types for biological network plot
#
#    Arguments:     None
#    Returned value:    None
#
#    Example:   supportedBioNetCircosPlotType()
#
#    Last revised on July 14, 2014
#

showSupportedBioNetCircosPlotType <- function() {

    cat("Currently supported plot type:\n\n")
    cat("1. group (for sample inforamtion, such as case/control),\n")
    cat("2. bar (for data values between 0 an 1),\n")
    cat("3. points (such as copy number variation),\n")
    cat("4. heatmap (such as gene expression),\n")
    cat("5. line\n")
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Create color map for heatmap plot. This is adopted from RCircos package
#
#    Arguments:
#
#        colorType: character vector, one of following:
#
#            BlueWhiteRed:  colors from blue to white then red
#            GreenWhiteRed: colors from green to white then red
#            GreenYellowRed: colors from green to yellow then red
#            GreenBlackRed: colors from green to black then red
#            YellowToRed:   colors from yellow to red
#            BlackOnly: default black colors
#
#    Example:   internal use only
#
#    Last revised on October 28, 2013 
#

getHeatmapColorScales<-function(colorType) {

    allOnes  <- seq(1, 1, length=256)
    allZeros <- seq(0, 0, length=256)
    one2zeor <- seq(1, 0, length=256)
    zero2one <- seq(0, 1, length=256)

    #   Blue, White, and Red
    #   +++++++++++++++++++++++++++++++++++++++++
    if(colorType == "BlueWhiteRed") {

        blueRamp <- rgb(zero2one, zero2one, allOnes)
        redRamp  <- rgb(allOnes, one2zeor, one2zeor)
        colorRamp <- cbind(blueRamp, redRamp)

    #   Green, White, and Red
    #   +++++++++++++++++++++++++++++++++++++++++
    } else if (colorType == "GreenWhiteRed") {

        redRamp   <- rgb(allOnes, one2zeor, one2zeor)
        greenRamp <- rgb(zero2one, allOnes, zero2one)
        colorRamp <- cbind(greenRamp, redRamp)

    #   Green, Yellow, and Red
    #   +++++++++++++++++++++++++++++++++++++++++
    } else if (colorType == "GreenYellowRed"){

        redRamp <- rgb(allOnes, one2zeor, allZeros)
        greenRamp <- rgb(zero2one, allOnes, allZeros)
        colorRamp <- cbind(greenRamp, redRamp)

    #   Green, Black, and Red
    #   +++++++++++++++++++++++++++++++++++++++++
    } else if (colorType == "GreenBlackRed"){

        redRamp <- rgb( zero2one, allZeros, allZeros)
        greenRamp <- rgb(allZeros, one2zeor, allZeros)
        colorRamp <- cbind(greenRamp, redRamp)

    #   Yellow to Red
    #   +++++++++++++++++++++++++++++++++++++++++
    } else if (colorType == "YellowToRed") {
    
        colorRamp <- rgb(allOnes, one2zeor, allZeros)

    #   black only
    #   +++++++++++++++++++++++++++++++++++++++++
    } else {
        colorRamp <- rgb(one2zeor, one2zeor, one2zeor)
        cat(paste("Warning: an unsupported color type",
                "defined and black will be used!\n"))
    }

    return (colorRamp)
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Plot a color scale for heatmap. This function is called by legend plot 
#    function and is not intended for direct call by users.
#
#    Arguments:
#
#        coorX:  numeric, x coordinates for the top left of color scale
#        coorY:  numeric, y coordinates for the top left of color scale
#
#        scaleWidth:  non-negative numeric, width of color scale
#        scaleHeight: non-negative numeric, height of color scale
#
#        colorType: character vector, one of following:
#
#            BlueWhiteRed:  colors from blue to white then red
#            GreenWhiteRed:  colors from green to white then red
#            GreenYellowRed: colors from green to yellow then red
#            GreenBlackRed:  colors from green to black then red
#            YellowToRed:    colors from yellow to red
#            BlackOnly:      default black colors
#
#        minValue:  The smallest value associated with the lowest color
#        maxValue:  The highest value associated with the highest color
#        direction: character, either "h" for horizotal or "v" for vertical
#
#    Returned value: None
#
#    Example:   plotHeatmapColorScale(coorX=100, coorY=-300, 
#                    scaleWidth=100, scaleHeight=10, 
#                    colorType="BlueWhiteRed", 
#                    minValue=-10, maxValue=10)
#
#   Last revised on September 24, 2014
#

plotHeatmapColorScale <- function(coorX, coorY, colorType="BlueWhiteRed", 
        scaleWidth, scaleHeight, minValue, maxValue, direction="h") {

    if(length(coorX) == 0 || length(coorY) == 0)
        stop("x and y coordinatges for color scale must be defined.")

    if(length(scaleWidth) == 0 || length(scaleHeight) == 0 || 
        scaleWidth<0 || scaleHeight<0)
        stop("scaleWidth and scaleHeight must be defined as non-negative.")

    if(is.null(minValue) || is.null(maxValue) ) {
        stop("minValue and maxValue must be defined.")
    }

    direction <- tolower(direction)
    if(direction != "h" && direction != "v") 
        stop("Color scale direction should be h or v")

    colorRamp <- getHeatmapColorScales(colorType)
    totalRect <- nrow(colorRamp)*ncol(colorRamp)

    if(direction == "h") {

        rectWidth <- scaleWidth/totalRect
        rectHeight <- scaleHeight

        yTop <- coorY
        yBottom <- yTop - scaleHeight

        for(aRect in 1:totalRect) {

            xLeft <-  coorX + (aRect-1)*rectWidth
            xRight <- xLeft + rectWidth

            rect(xLeft, yBottom, xRight, yTop, col=colorRamp[aRect],  
                        border = NA)
        }

        text(coorX, coorY-(scaleHeight/2), minValue, pos=2)
        text(coorX+scaleWidth, coorY-(scaleHeight/2), maxValue, pos=4)

    } else {
        rectWidth  <- scaleHeight
        rectHeight <- scaleWidth/totalRect

        xLeft <-  coorX
        xRight <- xLeft + rectWidth

        for(aRect in totalRect:1) {

            yTop <- coorY - (aRect-1)*rectHeight
            yBottom <- yTop - rectHeight
            rect(xLeft, yBottom, xRight, yTop, col=colorRamp[aRect],  
                border = NA)
        }

        text(coorX+(rectWidth/2), coorY, minValue, pos=3)
        text(coorX+(rectWidth/2), coorY-scaleHeight, maxValue, pos=1)
    }
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Convert numeric matrix to z-scores by row
#
#    Argument: a data frame with first column as row ID
#    Returned value:    A data frame with z scores for each row
#
#    Example:        zscores <- getZScores(exprData);
#
#    Last revised on Dec 16, 2014
#

convertToZScores <- function(exprData) {

    if(is.data.frame(exprData) == FALSE) 
        stop("Input data must be data frame with first column as row ID.")

    zscore <- as.matrix(exprData[,-1])
    for(aRow in 1:nrow(exprData)) {
        theExpr <- as.numeric(zscore[aRow,]);
        aMean <- mean(theExpr);
        aSD   <- sd(theExpr);
        zscore[aRow,] <- (theExpr-aMean)/aSD;
    }

    return (data.frame(Gene=as.character(exprData[,1]), zscore));
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Add legend to biomatrix layout to show heatmap scale, category definition,
#    and binary data definition.  
#
#    Arguments:
#
#        heatmapNames:  character vector, names of heatmap data
#        categoryNames: character vector, definition of categories
#        bionaryNames:  character vector, definition of binary items
#        HeatmapMin:    numeric, minimum value for heatmap
#        HeatmapMax:    numeric, maximum value for heatmap
#        colorType: character vector, default: "BlueWhiteRed", 
#                    other valid values are "GreenWhiteRed", "GreenYellowRed",
#                    "GreenBlackRed", "YellowToRed", and "BlackOnly"
#
#    Returned value: none
#    Example: bioMatrixLegend(heatmapNames=c("RNA", "miRNA"),
#            categoryNames=c("Methylation High", "Methylation Low")
#            binaryNames=c("DNA Amplification", "DNA deletion") )
#
#
#   Last revised on: January 30, 2015
#
#

bioMatrixLegend <- function(heatmapNames=NULL, categoryNames=NULL, 
            binaryNames=NULL, heatmapMin=-3, heatmapMax=3,
            colorType="BlueWhiteRed") {

    if(length(heatmapNames)>0) {

        colScalePosX <- getBioMatrixGeneLabelWidth();
        colScalePosY <- 0.25;
        scaleLength <- floor(getBioMatrixDataAreaWidth()/3);

        plotHeatmapColorScale(coorX=colScalePosX, coorY=colScalePosY, 
            colorType=colorType, scaleWidth=scaleLength, scaleHeight=0.4,
            minValue=heatmapMin, maxValue=heatmapMax);

            textPosX <- getBioMatrixGeneLabelWidth();
            textPosY <- -0.45;

            if (length(heatmapNames) == 1) {
                    scaleText <- heatmapNames;
            } else if(length(heatmapNames) == 2) {
                    scaleText <- paste("Top: ", heatmapNames[1], "  ", 
                                    "Bottom: ", heatmapNames[2])
            } else { stop("Incorrect heatmap names defined.") }

            text(textPosX, textPosY, scaleText, pos=4, offset=0)
    }

    if(length(categoryNames)>0) {
        legendX <- floor(getBioMatrixDataAreaWidth()/2) + 1
        legendY <- 0.5
        borderColors <- rev(getCaOmicsVColors())
        legend(legendX, legendY, legend=categoryNames, fill="white",  
            border=borderColors[1:length(categoryNames)], bty="n",
            xjust=0, yjust=1)
    }


    if(length(binaryNames)>0) {
        legendX <- floor(getBioMatrixDataAreaWidth()/4*3) + 2
        legendY <- 0.5
        caOmicsVColors <- getCaOmicsVColors()
        legend(legendX, legendY, legend=binaryNames, pch=19, 
            col=c(caOmicsVColors[4], caOmicsVColors[3]), 
            bty="n", xjust=0, yjust=1)
    }
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Add legend to bioNetCircos layout to show heatmap scale, category 
#    definition, and binary data definition. 
#
#    Arguments:
#
#        dataNames:     character vector, names of dataset(s) to be ploted
#        textCoor:      numeric vector of length 2, x and y coordinates to
#                       plot text
#        heatmapCoor:   numeric vector of length 2, x and y coordinates to
#                      plot color scale
#        scaleWidth:    non-negative numeric , width of color scale in inch
#        scaleHeight:   non-negative numeric , height of color scale in inch 
#        minHeatmap:    numeric, minimum value for heatmap, default -3
#        maxHeatmap:    numeric, maximum value for heatmap, default 3 
#        colorType:     character vector, default: "BlueWhiteRed", other valid 
#                       values are "GreenWhiteRed", "GreenYellowRed",
#                       "GreenBlackRed", "YellowToRed", and "BlackOnly"
#        direction:     character, either "h" for horizontal or "v" for vertical
#
#   Returned value: none
#
#   Example: bioMatrixLegend(dataNames=c("RNASeq", "miRNASeq",
#                    "Methylation", "DNA Amplification"), textCoor=c(10, 10)
#                    scaleWidth=4, scaleHeight=0.25)
#
#   Last revised on: January 30, 2015
#
#

bioNetLegend <- function(dataNames, textCoor=NULL, heatmapCoor=NULL, 
            scaleWidth, scaleHeight, heatmapMin=-3, heatmapMax=3, 
            colorType="BlueWhiteRed", direction="h") {

    bioNetGraph <- getBioNetGraph()
    if(is.null(bioNetGraph$layout) == TRUE)
        stop("Node layout has not been initialized.\n") 

    RangeX <- range(bioNetGraph$layout[,1])
    RangeY <- range(bioNetGraph$layout[,2])

    if(length(heatmapCoor) != 2) {
        heatmapCoor <- c(RangeX[2], RangeY[2])
        scaleWidth <- RangeX[2]*0.25
        scaleHeight <- RangeY[2]*0.05
    }
    
    plotHeatmapColorScale(coorX=heatmapCoor[1], 
            coorY=heatmapCoor[2], colorType=colorType, 
            scaleWidth=scaleWidth, scaleHeight=scaleHeight, 
            minValue=heatmapMin, maxValue=heatmapMax, 
            direction=direction);

    if(length(textCoor) !=2 ) 
        textCoor <- c(RangeX[2], RangeY[2]-scaleHeight*2)

    text(textCoor[1], textCoor[2], "From center to outer:", pos=4)
    legend(textCoor[1], textCoor[2], legend=dataNames, pch=19, 
            bty="n", xjust=0, yjust=1);
}

#   Last modified on May 18, 2015
#   __________________________________________________________________________
