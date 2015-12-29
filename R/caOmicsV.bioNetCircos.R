#
#    bioNetCircos layout plot 
#
#    Data types to  be visualized included:
#
#    Clinical data (tissue type, diagnosis, metastasis, ...
#    Gene expression (up- and down-regulations from microarray or NGS)
#    Copy number variations (deletions, insertions, and amplifications)
#    Mutations
# 
#    Date created: June 17, 2014
#    Revised on May 18, 2015 in compliance with Bioconductor coding style
# 
#    Hongen Zhang, Ph.D. (hzhang@mail.nih.gov)
#
#    Genetics Branch
#    Center for Cancer Research 
#    National Cancer Institute
#    National Institutes of Health
#    Bethesda, Maryland 20892
#





#    __________________________________________________________________________
#    **************************************************************************
#
#    Initialize the plot parameters, default coordinates of a circular line,  
#    and node layout for BioNetCircos plot.  
#
#    Arguments:
#
#        bioNet:        An igraph object representing a biological network or
#                       a customized igraph object
#        totalSamples:  non-negative integer, total number of samples, must
#                       defined by user
#        sampleWidth:   non-negative integer, Points per sample along the node
#                       circumference, default 100
#        nodeRadius:    non-negative numeric, Radius of the node, default 1
#        nodePadding:   non-negative numeric, minimum space between plot area
#                       of two nodes, relative to node radius, default 1
#        plotAreaWidth: non-negative numeric, outside boundary of plot area of
#                       a node, relative to node radius, default 1
#        layout:        layout of igraph object
#
#
#    Returned value:    None. Key parameters are stored in CA_OMICS_ENV space
#
#    Example:   initializeBioNetCircos(bioNet)
#             initializeBioNetCircos(bioNet, totalSamples=100, 
#                        sampleWidth=100)
#
#    Last revisited on August 12, 2014
#
#

initializeBioNetCircos <- function(bioNet, totalSamples=100, sampleWidth=100,
            nodeRadius=1, nodePadding=1, plotAreaWidth=1, 
            layout=layout.fruchterman.reingold(bioNet)) {

    if (!is.igraph(bioNet)) 
        stop("The first argument must be an igraph object!")

    if(totalSamples<0 || sampleWidth<0 || nodeRadius<0 ||
            nodePadding<0 || plotAreaWidth<0 ) {
        stop("Negative value is not allowed for argument(s)") 
    }

    if((sampleWidth %% 2) == 1)  sampleWidth<-sampleWidth + 1

    setBioNetPlotParameters(totalSamples, sampleWidth, 
            nodeRadius, nodePadding, plotAreaWidth)

    setBioNetCircosBasePlotPositions(totalSamples, sampleWidth)
    setBioNetNodeLayout(bioNet, layout)
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    Get x and y coordinates for points on a circular line. These coordinates
#    are relative to the point (0, 0) and will be transformed for different 
#    nodes. The default radius of the circle is 1.
#
#    Arguments:
#
#        totalSamples: non-negative integer, total number of samples
#        sampleWidth:  non-negative integer, points of a sample on the
#                      circumference of a node 
#
#    Returned value: None
#    Example:        setBioNetCircosBasePlotPositions(155, 100)
#
#   Last revisited on June 24, 2014
#

setBioNetCircosBasePlotPositions <- function(totalSamples=100, sampleWidth=100)
{
    if(totalSamples<0 || sampleWidth<0) 
            stop("Arguments must be non-negative integer.\n")

    totalPoints <- totalSamples * sampleWidth
    if(totalPoints<10000) {  
        sampleWidth <- ceiling(10000/totalSamples)
        if((sampleWidth %% 2) == 1) sampleWidth<-sampleWidth + 1 
        totalPoints <- totalSamples * sampleWidth
    }

    interval <- 2*pi/totalPoints
    baseVal <- seq(0, 2*pi, interval)

    corX <- sin(baseVal)
    corY <- cos(baseVal)

    degree <- rep(0, length(baseVal))
    mid <- round((length(baseVal)-1)/2, digits=0) + 1

    totalPoints <- length(baseVal)
    degree[1:mid] <- 90 - (baseVal[1:mid]*180/pi)
    degree[(mid+1):totalPoints] <- 270 - (baseVal[(mid+1):totalPoints]*180/pi)

    basePostions <- data.frame(corX, corY, degree)

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())
    caOmicsVEnvironment[["BioNet_Base_Location"]] <- basePostions
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Set up plot parameters for caOmicsV bioNetCircos layout.
#
#    Arguments:
#
#        totalSamples:  non-negative integer, total number of samples
#        sampleWidth:   non-negative integer, points of a sample on the
#                       circumference of a node 
#        nodeRadius:    non-negative numeric, radius of node on igraph
#        nodePadding:   non-negative numeric, empty area between two nodes
#        plotAreaWidth: non-negative numeric, outside boundary of plot area
#
#    Returned value:    none
#
#    Example:   setBioNetPlotParameters(totalSamples=100, sampleWidth=100, 
#                    nodeRadius=1, nodePadding=0.2, plotAreaWidth=10)
#
#    Last revisited on July 17, 2014
#
#

setBioNetPlotParameters <- function(totalSamples, sampleWidth, 
            nodeRadius, nodePadding, plotAreaWidth)
{
    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    caOmicsVEnvironment[["BioNet_Total_Sample"]]   <- totalSamples
    caOmicsVEnvironment[["BioNet_Sample_Width"]]   <- sampleWidth
    caOmicsVEnvironment[["BioNet_Node_Radius"]]    <- nodeRadius
    caOmicsVEnvironment[["BioNet_Node_Padding"]]   <- nodePadding
    caOmicsVEnvironment[["BioNet_PlotArea_Width"]] <- plotAreaWidth
    caOmicsVEnvironment[["BioNet_PlotArea_Inner"]] <- nodeRadius
    caOmicsVEnvironment[["BioNet_PlotArea_Outer"]] <- nodeRadius

    setCaOmicsVColors();
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Set up layout of nodes on the biological network. The layout is taken 
#    from a igraph layout and make necessary scaling to alocate circos plot 
#    area for each node
#
#    Prerequisite:  igraph package muse be loaded first
#
#    Arguments:
#            bioNet:    an igraph object representing a biological network
#            layout:    a two dimensional matrix of x and y coordinates for 
#                       node centers
#
#    Returned value:    None. 
#
#    Example:   layout <- layout.fruchterman.reingold.grid(bioNet)
#               setBioNetNodeLayout(bionet, layout)
#
#    Last revisited on July 29, 2014
#
#

setBioNetNodeLayout <- function(bioNet, layout=layout.auto(bioNet)) {

    if (!is.igraph(bioNet)) 
        stop("An igraph object is required for layout setting.") 

    nodeRadius  <- getBioNetNodeRadius()
    nodePadding <- getBioNetNodePaddingScale()
    plotAreaWid <- getBioNetPlotAreaWidth()

    layout   <- layout.norm(layout, -1, 1, -1, 1)
    minDist  <- min(dist(layout))
    minSpace <- nodeRadius*(2 + plotAreaWid*2 + nodePadding)
    layout   <- layout*(minSpace/minDist) 

    bioNet <- set.graph.attribute(bioNet, "layout", layout)

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())
    caOmicsVEnvironment[["BioNet_Graph"]] <- bioNet
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Record the plotted area boundary for all nodes on the igraph. These 
#    boundary may be needed for drawn customized arrows and label node names.
#
#    Argument:
#
#        inner: non-negative numeric, the inner boundary of area that has been 
#               plotted
#        outer: non-negative numeric, the outer boundary of area that has been 
#               plotted
#
#    Return value:  None.
#
#    Example:   resetBioNetNodePlotAreaBoundary(1, 1.5)
#
#    Last revised on August 19, 2014
#
#

resetBioNetNodePlotAreaBoundary <- function(inner=getBioNetNodeRadius(), outer)
{
    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    lastInner <- caOmicsVEnvironment[["BioNet_PlotArea_Inner"]]
    if(inner<lastInner)
        caOmicsVEnvironment[["BioNet_PlotArea_Inner"]] <- inner

    lastOuter <- caOmicsVEnvironment[["BioNet_PlotArea_Outer"]]
    if(lastOuter<outer)
        caOmicsVEnvironment[["BioNet_PlotArea_Outer"]] <- outer 
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Get methods to retrieving caOmicsV objects stored in caOmicsV environment
#
#    Argument:     None
#    Return value: one of objects stored in caOmicsV environment.
#
#    Example:   plotPositions <- getBasePositions()
#
#    Last revised on AUgust 12, 2014
#

getBioNetPlotTotalSample <-function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioNet_Total_Sample"]])
}

getBioNetPlotSampleWidth<-function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioNet_Sample_Width"]])
}

getBioNetNodeRadius <-function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioNet_Node_Radius"]])
}

getBioNetNodePaddingScale <-function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioNet_Node_Padding"]])
}

getBioNetPlotAreaWidth <-function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioNet_PlotArea_Width"]])
}

getBioNetBasePositions<-function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioNet_Base_Location"]])
}

getBioNetGraph<-function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    return (caOmicsVEnvironment[["BioNet_Graph"]])
}

getBioNetNodePlotAreaBoundary <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    inner <- caOmicsVEnvironment[["BioNet_PlotArea_Inner"]]
    outer <- caOmicsVEnvironment[["BioNet_PlotArea_Outer"]]

    return (c(inner=inner, outer=outer))
}

getBioNetNodeParameters <- function() {

    caOmicsVEnvironment <- NULL
    caOmicsVEnvironment <- get(CA_OMICS_NAME, envir=globalenv())

    totalSamples  <- caOmicsVEnvironment[["BioNet_Total_Sample"]]
    sampleWidth   <- caOmicsVEnvironment[["BioNet_Sample_Width"]]
    nodeRadius    <- caOmicsVEnvironment[["BioNet_Node_Radius"]]
    nodePadding   <- caOmicsVEnvironment[["BioNet_Node_Padding"]]
    plotAreaWidth <- caOmicsVEnvironment[["BioNet_PlotArea_Width"]]

    inner <- caOmicsVEnvironment[["BioNet_PlotArea_Inner"]]
    outer <- caOmicsVEnvironment[["BioNet_PlotArea_Outer"]]

    return (bioNetParams=list(totalSamples, sampleWidth, nodeRadius, 
            nodePadding, plotAreaWidth, inner, outer))
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    Display biological network before plot omics data on each node. 
#
#    Vertex will be plotted as colored polygons
#    Vertex size is controlled by node radius 
#
#    Prerequisite:  the igraph object must have a layout attached
#
#    Arguments:    character vector for color or a R color specification
#    Return value: None
#
#    Example:        showBioNetNodesLayout(bgColor=grey(0.75, alpha=0.5))
#
#    Last revisited on July 23, 2014
#
#

showBioNetNodesLayout <- function(bgColor=grey(0.75, alpha=0.5)) {

    bioNetGraph <- getBioNetGraph()
    nodeRadius  <- getBioNetNodeRadius()
    plotPositions <- getBioNetBasePositions()

    if(is.null(bioNetGraph$layout) == TRUE)
        stop("Node layout has not been initialized.\n") 

    RangeX <- range(bioNetGraph$layout[,1])
    RangeY <- range(bioNetGraph$layout[,2])

    RangeX[2] <- RangeX[2]*1.5

    plot(bioNetGraph, rescale=FALSE, vertex.label=NA, xlim=RangeX, ylim=RangeY)

    text(bioNetGraph$layout[,1], bioNetGraph$layout[,2], V(bioNetGraph))

    setBioNetPlotAreaBackground(bgColor)
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Change the plot area background of igraph node. Use white color to erase 
#    background and and grey to show the plot area boundary. This will not 
#    change the node layout
#
#    Argument:     character vector for R color or a R color specification
#    Return value: None
#
#    Example:   setBioNetPlotAreaBackground(bgColor=red(0.75, alpha=0.5))
#
#    Last updated on August 19, 2014
#

setBioNetPlotAreaBackground <- function(bgColor=grey(0.75, alpha=0.5)) {

    bioNetGraph <- getBioNetGraph()
    if(is.null(bioNetGraph$layout) == TRUE)
        stop("Node layout has not been initialized.\n") 

    nodeRadius  <- getBioNetNodeRadius()
    plotAreaWid <- getBioNetPlotAreaWidth()
    plotRadius  <- nodeRadius*(1 + plotAreaWid)

    plotPositions <- getBioNetBasePositions()
    inPositions   <- plotPositions[,1:2]*nodeRadius
    outPositions  <- plotPositions[,1:2]*plotRadius

    end <- nrow(inPositions)
    polygonX <- c(outPositions[1:end,1], inPositions[end:1,1])
    polygonY <- c(outPositions[1:end,2], inPositions[end:1,2])

    rgb.val <- as.vector(col2rgb(bgColor))/255
    bgColor <- rgb(red=rgb.val[1], green=rgb.val[2], 
                    blue=rgb.val[3], alpha=0.5)

    for(nodeIndex in seq_len(nrow(bioNetGraph$layout)))  {

        nodeCenter <- bioNetGraph$layout[nodeIndex,]
        areaX <- polygonX + nodeCenter[1]
        areaY <- polygonY + nodeCenter[2]

        polygon(areaX, areaY, col=bgColor, border=NA)
    }
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Erase the node background for all nodes, usually after node names are 
#    plotted and before plotting data tracks. The area to be erased has same
#    radius as the node
#
#    Arguments:     None
#    Return value:  None
#
#    Example:   eraseBioNetNode()
#
#    Last revisited on July 23, 2014
#

eraseBioNetNode <- function() {

    bioNetGraph <- getBioNetGraph()
    nodeRadius  <- getBioNetNodeRadius()
    plotAreaWid <- getBioNetPlotAreaWidth()
    plotRadius  <- nodeRadius*(1 + plotAreaWid)

    plotPositions <- getBioNetBasePositions()
    nodeArea <- plotPositions[,1:2]*plotRadius

    for(nodeIndex in seq_len(nrow(bioNetGraph$layout)))  {

        nodeCenter <- bioNetGraph$layout[nodeIndex,]
        nodeX <- nodeArea[,1] + nodeCenter[1]
        nodeY <- nodeArea[,2] + nodeCenter[2]

        polygon(nodeX, nodeY, col="white", border="white")
    }
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    Calculate x and y coordinates for sample plot positions on default node.  
#    The output will be a three column matrix representing the left, center, 
#    and right position for each sample on circumference of default node. The 
#    center positions are for points plot and others are for polygon plot. 
#    This function is for internal use.
#
#    Arguments:    totalSamples, non-negative integer, total number of samples
#    Return value: matrix with index of x and y coordinates for each sample. 
#
#    Example:   pointIndex <- getPlotCoordinates(totalSamples=100)
#
#    Last revisited on July 7, 2014
#

getBioNetSamplePlotPosition <- function(totalSamples) {

    if(is.numeric(totalSamples) == FALSE || totalSamples<0 )
            stop("Incorrect total samples defined!\n") 

    basePositions <- getBioNetBasePositions()
    totalPoints   <- nrow(basePositions) - 1

    interval <- totalPoints/totalSamples
    shiftBy  <- interval/2

    right    <- interval*c(1:totalSamples)+1
    center   <- right - shiftBy
    left     <- center - shiftBy
    sampleLocation <- cbind(left, center, right)

    return (sampleLocation)
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    Calculate plot positions including node center, outer boundary, inner 
#    boundary, and index of point coordinates.
#
#    Arguments:
#
#    nodeCenter: x and y coordinates of the node center
#    outer:      non-negative numeric, outer limit of plot track and relative
#                node center (0,0)
#    inner:      non-negative numeric, inner limit of plot track and relative
#                to node center (0,0)
#
#    All arguments should be validated in advance.
#
#    Returned value: a list that contains node centers, outer boundary, inner 
#                    boundary, and index of point coordinates for polygons.
#
#    Example: locations <- getBioNetPlotLocations(nodeCenter=5,  
#                                    outer=1.9, inner=1.5)
#
#   Last revised on:    February 2, 2015
#
#

getBioNetPlotLocations <- function(nodeCenter, outer, inner) {

    theCenter <- as.numeric(nodeCenter)
    theRadius <- getBioNetNodeRadius()

    basePositions <- getBioNetBasePositions()
    outLocations  <- basePositions[,1:2]*theRadius*outer
    inLocations   <- basePositions[,1:2]*theRadius*inner 

    totalSamples <- getBioNetPlotTotalSample()
    theIndex <- getBioNetSamplePlotPosition(totalSamples)

    return (list(nodeCenter=theCenter, outPositions=outLocations, 
            inPositions=inLocations, positionIndex=theIndex))
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    The main plot function for caOmicsV bioNetCircos layout. It plots one  
#    type data on all nodes of the network.
#
#    Arguments:
#
#    dataValues:    nemeric matrix, the data used for generatrion of  
#                    biological network
#    plotType:      character vector, must be one of "group", "bar",  
#                    "scatters", "heatmap", "line"
#    outer:         non-negative numeric, the boundary of plot area far 
#                    from node center
#    inner:         non-negative numeric, the boundary of plot area close 
#                    to node center
#    plotColors:    character vector or vectors fo R color specification, 
#                    colors for each sample
#    maxValue:       numeric, the biggest value of plot data
#    minValue:       numeric, the smallest  value of plot data
#
#    Return value:  None
#
#    Example:      plotColors <- c(rep("red", 20), rep("blue", 30))
#                bioNetCircosPlot(dataValues=p53, plotType="polygon", 
#                        outer=0.9, inner=0.7, plotColors=plotColors)
#
#    Last revisited on August 14, 2014
#
#

bioNetCircosPlot <- function(dataValues=NULL, plotType="polygon", outer,  
            inner, plotColors=NULL, maxValue=NULL, minValue=NULL) {

    supportedType <- getCaOmicsVPlotTypes();
    plotType <- tolower(plotType);

    
    if(!plotType %in% supportedType)  stop("Unsupported plot type!\n")  
    if(is.null(dataValues))  stop("Plot data is missing.")  

    if(is.vector(plotColors) == FALSE || is.character(plotColors) == FALSE)
            stop("Plot colors must be held with vector!\n") 
    if(is.null(dataValues))  stop("Plot data is missing.")  

    if(is.numeric(outer) == FALSE || is.numeric(inner) == FALSE) 
        stop("Plot area boundary must be numeric!\n")  
    if(outer<0 || inner<0) stop("Plot area boundary cannot be negative.\n")  

    if(outer == inner) stop("Plot area is too small.\n") 
    if(outer<inner) { temp <- outer; outer <- inner; inner <- temp }

    bioNetGraph <- getBioNetGraph()
    if(is.null(bioNetGraph$layout)) 
        stop("Node layout has not been initialized.\n") 


    #   plot data as polygons, no color definition need
    #   ++++++++++++++++++++++++++++++++++++++++++++++++++

    if(plotType == supportedType[1]) {
        plotBioNetPolygons(dataValues, outer, inner) 

    #   plot data as bars
    #   ++++++++++++++++++++++++++++++++++++++++++++++++++

    } else if(plotType == supportedType[2]) {
        plotBioNetBars(dataValues, outer, inner, plotColors)

    #   plot data as points
    #   ++++++++++++++++++++++++++++++++++++++++++++++++++

    } else if(plotType == supportedType[3]) {
        plotBioNetPoints(dataValues, outer=outer, inner=inner)

    #   plot data as heatmap
    #   ++++++++++++++++++++++++++++++++++++++++++++++++++

    } else if(plotType == supportedType[4]) {
        plotBioNetHeatmap(dataValues, maxValue, minValue, 
                    outer, inner, plotColors)

    #   plot data as lines between two neighbor data points
    #   ++++++++++++++++++++++++++++++++++++++++++++++++++

    } else if(plotType == supportedType[5]) {
        plotBioNetLines(dataValues, maxValue, minValue, 
                    outer, inner, plotColors)

    #   plot category data as colored polygons
    #   ++++++++++++++++++++++++++++++++++++++++++++++++++

    } else if(plotType == supportedType[6] || plotType == supportedType[7]) {
        plotBioNetPolygons(dataValues, outer, inner)  

    #   error report
    #   ++++++++++++++++++++++++++++++++++++++++++++++++++

    } else { stop("Nothing could be done!\n") }

    resetBioNetNodePlotAreaBoundary(inner, outer)
}




#    ________________________________________________________________________
#    ************************************************************************
#
#    Plot group data such as clinical Features on each node. By passing a  
#    matrix with more than one row, this function could be used to plot 
#    category/binary data.
#
#    Arguments:
#
#        dataValues:    matrix of character or numeric for category data
#        outer:         non-negative numeric, the outer boundary of plot area
#        inner:         non-negative numeric, the inner boundary of plot area
#
#    Return value:  None
#
#    Example:   plotColors <- c(rep("red", 20), rep("blue", 30))
#               plotBioNetGroupData(dataValues, outer=0.9, inner=0.7)
#
#    Last revisited on August 14, 2014
#

plotBioNetPolygons <- function(dataValues, outer, inner) {

    bioNetGraph <-  getBioNetGraph()

    groupNames <- unique(as.vector(dataValues));
    numOfGroup <- length(groupNames);

    colorSet <- getCaOmicsVColors()
    if(numOfGroup<=length(colorSet)) {
        groupColors <- colorSet[1:numOfGroup] 
    } else { 
        groupColors <- rainbow(numOfGroup)  
    }

    for(nodeIndex in seq_len(nrow(bioNetGraph$layout)))  {

        nodeCenter <- as.numeric(bioNetGraph$layout[nodeIndex,])
        plotLocations <- getBioNetPlotLocations(nodeCenter, outer, inner)

        if(nrow(dataValues) == 1)  { 
            rowIndex <- 1
        } else { rowIndex <- nodeIndex }

        plotData <- dataValues[rowIndex,]
        plotColors <- rep(groupColors[1], length(plotData));

        for(aGroup in 1:numOfGroup) {  
            items <- which(plotData == groupNames[aGroup])
            plotColors[items] <-  groupColors[aGroup] 
        }

        totalSamples <- getBioNetPlotTotalSample()
        for(a.sample in seq_len(totalSamples)) {
            start <- plotLocations$positionIndex[a.sample, 1]
            end   <- plotLocations$positionIndex[a.sample, 3]
            theColor <- plotColors[a.sample]

            polygonX <- c(plotLocations$outPositions[start:end,1], 
                            plotLocations$inPositions[end:start,1])
            polygonY <- c(plotLocations$outPositions[start:end,2], 
                            plotLocations$inPositions[end:start,2])

            polygonX <- polygonX + nodeCenter[1]
            polygonY <- polygonY + nodeCenter[2]

            polygon(polygonX, polygonY, col=theColor, border=NA)
        }
    }
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    Bar plot for each sample, mostly used for displaying percentages such as 
#    coverage from NGS data or ratios of expressed genes of a pathway ... 
#
#    Arguments:
#
#    dataValues: numeric matrix with range of 0 ~ 1 for bar height. total 
#                rows of the matrix must be same as the number of nodes 
#                and row names must be same as the vertex names in bioNetGraph
#    outer:      non-negative numeric, the outer boundary of plot area 
#    inner:      non-negative numeric, the inner boundary of plot area
#    plotColors: colors for each sample
#
#    Return value:  None
#
#    Example:   plotBioNetFrequency(dataValues, nodeIndex, outer, inner, 
#                                  plotColors)
#
#    Last revisited on August 14, 2014
#

plotBioNetBars <- function(dataValues, outer, inner, plotColors) {

    if(ncol(dataValues) !=length(plotColors)) { 

        cat("length of data values and colors are different.\n")
        cat("First color will be used for all samples.\n")

        plotColors <- rep(plotColors[1], length(dataValues))  
    }

    barData <- as.vector(dataValues);
    if(is.numeric(barData) == FALSE) 
        stop("Bar plot data must be numeric between 0 and 1.")

    if(max(barData)>1 || min(barData)<0) {

        maxValue <- max(barData);
        minValue <- min(barData);
        dataRange <- maxValue - minValue;
        for(aRow in seq_len(nrow(dataValues)))
            dataValues[aRow,] <- (dataValues[aRow,]- minValue)/dataRange
    }

    bioNetGraph <-  getBioNetGraph()
    trackHeight <- outer-inner

    vertexNames  <- V(bioNetGraph)$name
    dataRowNames <- rownames(dataValues)

    for(nodeIndex in seq_len(nrow(bioNetGraph$layout))) {

        if(length(dataRowNames) == 1) dataIndex <- 1
        else {
            dataIndex <- grep(vertexNames[nodeIndex], dataRowNames)
            if(length(dataIndex) == 0) next
        }

        plotData   <- as.numeric(dataValues[dataIndex,])
        nodeCenter <- as.numeric(bioNetGraph$layout[nodeIndex,])
        plotLocations <- getBioNetPlotLocations(nodeCenter, outer, inner)

        for(a.sample in seq_len(length(plotData))) {

            if(plotData[a.sample] == 0) next 

            barHeight <- (trackHeight*plotData[a.sample])/inner
            theColor <- plotColors[a.sample]

            start <- plotLocations$positionIndex[a.sample, 1]
            end   <- plotLocations$positionIndex[a.sample, 3]

            polygonX <- c(plotLocations$inPositions[start:end,1]*(1+barHeight),
                            plotLocations$inPositions[end:start,1])
            polygonY <- c(plotLocations$inPositions[start:end,2]*(1+barHeight),
                            plotLocations$inPositions[end:start,2])

            polygonX <- polygonX + nodeCenter[1]
            polygonY <- polygonY + nodeCenter[2]

            polygon(polygonX, polygonY, col=theColor, border=NA)
        }
    }
}






#    __________________________________________________________________________
#    **************************************************************************
#
#    Heatmap plot for all nodes of an igraph. 
#
#    Arguments:
#
#        dataValues:    numeric matrix of log2 values. total rows of the matrix
#                       must be same as the number of nodes and row names must
#                       be same as the vertex names in bioNetGraph
#        maxValue:      numeric, maximum value for highest color in heatmap, 
#                       set to NULL to use the maximum value in expression 
#                       dataset
#        minValue:      numeric, minimum value for lowest color in heatmap, set
#                       to NULL to use the minimum value in expression dataset
#        outer:         non-negative numeric, the outrt boundary of plot area
#        inner:         non-negative numeric, the innner boundary of plot area
#        plotColors:    color map, one of "BlueWhiteRed", "GreenWhiteRed",  
#                       "GreenYellowRed", "GreenBlackRed", "YellowToRed", 
#                       "BlackOnly".
#
#    Return value:  None
#
#    Example: plotBioNetHeatmap(dataValues=p53, nodeIndex=3, outer=1.4, 
#                            inner=1.2, plotColors="BlueWhiteRed")
#
#    Last revised on August 14, 2014
# 
#

plotBioNetHeatmap <- function(dataValues, maxValue=NULL, minValue=NULL,  
            outer, inner, plotColors) {

    errMSG <- "Heatmap colors should be one from:

        BlueWhiteRed, 
        GreenWhiteRed, 
        GreenYellowRed, 
        GreenBlackRed, 
        YellowToRed.

        Please redefine the plotColors."

    if(length(plotColors)>1)  stop(errMSG)

    if(is.null(minValue) == FALSE && is.null(maxValue) == FALSE)  
        dataRange <- c(minValue, maxValue)  
    else
        dataRange <- c(min(dataValues), max(dataValues))

    colorRamp <- getHeatmapColorScales(plotColors)
    colorLevel <- seq(dataRange[1], dataRange[2], length=length(colorRamp))

    bioNetGraph <-  getBioNetGraph()
    vertexNames  <- V(bioNetGraph)$name
    dataRowNames <- rownames(dataValues)

    for(nodeIndex in seq_len(nrow(bioNetGraph$layout))) {

        if(length(dataRowNames) == 1) dataIndex <- 1
        else {
            dataIndex <- grep(vertexNames[nodeIndex], dataRowNames)
            if(length(dataIndex) == 0) next
        }

        plotData  <- as.numeric(dataValues[dataIndex,])
        nodeCenter <- as.numeric(bioNetGraph$layout[nodeIndex,])
        plotLocations <- getBioNetPlotLocations(nodeCenter, outer, inner)

        sampleColors <- rep(colorRamp[1], length(plotData))
        for(a.sample in seq_len(length(plotData))) {

            if(is.na(plotData[a.sample])) {
                sampleColors[a.sample] <- "gray"; next 
            }

            the.level <- which(colorLevel>=plotData[a.sample])
            sampleColors[a.sample] <- colorRamp[min(the.level)]
        }

        for(a.sample in seq_len(length(plotData)))  {

            start <- plotLocations$positionIndex[a.sample, 1]
            end   <- plotLocations$positionIndex[a.sample, 3]
            theColor <- sampleColors[a.sample]

            polygonX <- c(plotLocations$outPositions[start:end,1], 
                            plotLocations$inPositions[end:start,1])
            polygonY <- c(plotLocations$outPositions[start:end,2], 
                            plotLocations$inPositions[end:start,2])

            polygonX <- polygonX + nodeCenter[1]
            polygonY <- polygonY + nodeCenter[2]

            polygon(polygonX, polygonY, col=theColor, border=NA)
        }
    }
}






#    __________________________________________________________________________
#    **************************************************************************
#
#    Point plot for all nodes of a igraph. Character type and size could be   
#    set from the console since points() function is used.
#
#    Arguments:
#
#        dataValues: numeric matrix of plot data
#        maxValue:   numeric, the biggest value of plot data
#        minValue:   numeric, the smallest  value of plot data
#        outer:      non-negative numeric, the outer boundary of plot area
#        inner:      non-negative numeric, the inner boundary of plot area
#        plotColors: character vector or R color specification for each sample
#        sizeByValue: logic, use data value for point size (cex)
#        pch: same as pch parameter in par()
#
#    Return value:  None
#
#    Example:   plotBioNetPoints(dataValues, 1.5, 2)
#
#    Last revised on AUgust 14, 2014
#

plotBioNetPoints <- function(dataValues, maxValue=NULL, minValue=NULL, 
        outer, inner, plotColors=rep("black", ncol(dataValues)), 
        sizeByValue=FALSE, pch=".")  {

    if(ncol(dataValues) != length(plotColors))
        stop("Length of data values and colors must be same.\n") 

    dataLow <- min(as.vector(dataValues))
    dataTop <- max(as.vector(dataValues)) 
    if(is.null(minValue) == FALSE) dataLow <- min(dataLow, minValue) 
    if(is.null(maxValue) == FALSE) dataTop <- max(dataTop, maxValue)

    plotHeight <- outer - inner
    dataHeight <- dataTop - dataLow
    pointHeights <-(dataValues-dataLow)/dataHeight*plotHeight/inner

    bioNetGraph <-  getBioNetGraph()
    vertexNames  <- V(bioNetGraph)$name
    dataRowNames <- rownames(dataValues)

    totalSamples  <- length(dataValues[1,])
    for(nodeIndex in seq_len(nrow(bioNetGraph$layout))) {

        if(length(dataRowNames) == 1) 
            dataIndex <- 1
        else 
            dataIndex <- nodeIndex

        nodeCenter <- as.numeric(bioNetGraph$layout[nodeIndex,])
        plotLocations <- getBioNetPlotLocations(nodeCenter, outer, inner)
        drawBioNetNodeBackground(plotLocations)

        pointIndex     <- plotLocations$positionIndex[,2]
        pointLocations <- plotLocations$inPositions[pointIndex,] 
        pointLocations <- pointLocations * (1+pointHeights[nodeIndex,])

        if(sizeByValue == TRUE) 
            pcex <-  pointHeights
        else  
            pcex <- rep(0.5, totalSamples)


        for(aSample in seq_len(totalSamples)) {

            theColor <- plotColors[aSample]

            pointX <- pointLocations[aSample, 1] + nodeCenter[1]
            pointY <- pointLocations[aSample, 2] + nodeCenter[2]
            points(pointX, pointY, col=theColor, pch=pch, cex=pcex[aSample])
        }
    }
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    Plot lines between two samples
# 
#    Arguments:
#
#        dataValues:    numeric matrix of plot data
#        outer:      non-negative numeric, the outer boundary of plot area
#        inner:      non-negative numeric, the inner boundary of plot area  
#        maxValue:   numeric, the biggest value of plot data
#        minValue:   numeric, the smallest  value of plot data  
#        plotColors: character vector or R color specification for each sample
#
#    Return value:  None
#
#    Example:   plotLines(bioNetGraph, dataValues, maxValue, minValue, 
#                       outer, inner, plotColors)
#
#    Last revised on July25, 2014
# 

plotBioNetLines <- function(dataValues, outer, inner, maxValue=NULL, 
            minValue=NULL, plotColors=rep("black", ncol(dataValues))) {

    if(ncol(dataValues) !=length(plotColors))
        stop("Length of data values and colors must be same.\n")

    totalSamples <- ncol(dataValues)

    dataLow <- min(dataValues)
    dataTop <- max(dataValues) 
    if(is.null(minValue) == FALSE) dataLow <- min(dataLow, minValue)  
    if(is.null(maxValue) == FALSE) dataTop <- max(dataTop, maxValue) 

    plotHeight <- outer - inner
    dataHeight <- dataTop - dataLow
    pointHeights <-(dataValues-dataLow)/dataHeight*plotHeight/inner

    bioNetGraph <-  getBioNetGraph()
    vertexNames  <- V(bioNetGraph)$name
    dataRowNames <- rownames(dataValues)

    for(nodeIndex in 1:nrow(bioNetGraph$layout)) {

        if(length(dataRowNames) == 1)  { 
                dataIndex <- 1
        } else {
            dataIndex <- grep(vertexNames[nodeIndex], dataRowNames)
            if(length(dataIndex) == 0) next
        }

        plotData  <- as.numeric(dataValues[dataIndex,])
        nodeCenter <- as.numeric(bioNetGraph$layout[nodeIndex,])

        plotLocations <- getBioNetPlotLocations(nodeCenter, outer, inner)
        drawBioNetNodeBackground(plotLocations)

        pointIndex     <- plotLocations$positionIndex[,2]
        pointLocations <- plotLocations$inPositions[pointIndex,] 
        pointLocations <- pointLocations* (1+pointHeights[nodeIndex,]) 

        for(a.sample in seq_len(totalSamples-1)) {

            theColor <- plotColors[a.sample]

            startX <- pointLocations[a.sample,1] + nodeCenter[1]
            endX <- pointLocations[a.sample+1,1] + nodeCenter[1]

            startY <- pointLocations[a.sample,2] + nodeCenter[2]
            endY <- pointLocations[a.sample+1,2] + nodeCenter[2]

            lines(c(startX, endX), c(startY, endY), col=theColor)
        }
    }
}






#    __________________________________________________________________________
#    **************************************************************************
#
#    Label node (vertex) names for customized locations. 
#
#    Argument:
#
#        nodeList:      no-negative numeric vector, which node(s) will be  
#                       labelled. Set NULL to label all nodes.
#        labelColor:    character vector or R color specification, colors for 
#                       node name(s)
#        lableLocation: character vector, location relative to node center,
#                       either "bottom", "left", "top", or "right"
#        labelOffset: non-negative numeric, distance from node outside boundary
#
#    Return value: None
#
#    Example:    labelBioNetNodeNames(1, "red", "bottom")
#                labelBioNetNodeNames(c(2:5), "blue", "right")
#
#    Last updated on August 20, 2014
#

labelBioNetNodeNames <- function(nodeList=NULL, labelColor="black",
    labelLocation=c("bottom", "left", "top", "right"), labelOffset=0.5) {

    nodeRadius  <- getBioNetNodeRadius()
    plotAreaWid <- getBioNetPlotAreaWidth()
    plotRadius  <- nodeRadius*plotAreaWid

    basePositions <- getBioNetBasePositions()
    outPositions  <- basePositions[,1:2]*plotRadius

    indexLength <- nrow(outPositions)
    labelLocation <- tolower(labelLocation)

    if(labelLocation == "bottom") {
        positionIndex <- round(indexLength/2, digits=0)
        nameLocation <- outPositions[positionIndex,]
        textPos <- 1
    } else if(labelLocation == "left") {
        positionIndex <- round(indexLength/4*3, digits=0)
        nameLocation <- outPositions[positionIndex,]
        textPos <- 2
    } else if (labelLocation == "top") {
        nameLocation <- outPositions[1,]
        textPos <- 3
    } else {
        positionIndex <- round(indexLength/4, digits=0)
        nameLocation <- outPositions[positionIndex,]
        textPos <- 4
    }

    bioNetGraph <- getBioNetGraph()
    nodeNames <- V(bioNetGraph)$name

    if(is.null(nodeList)) nodeList <- 1:nrow(bioNetGraph$layout) 

    for(nodeIndex in seq_len(length(nodeList))) {
        theNode <- nodeList[nodeIndex]
        nodeCenter <- bioNetGraph$layout[theNode,]

        nameX <- nameLocation[1] + nodeCenter[1]
        nameY <- nameLocation[2] + nodeCenter[2]
        text(nameX, nameY, nodeNames[theNode], pos=textPos, 
                offset=labelOffset, col=labelColor)
    }
}






#    __________________________________________________________________________
#    **************************************************************************
#
#    Draw background with grey or customized color for a plot track
#
#
#    Arguments:
#
#        trackLocations: an object returned from getBioNetPlotLocations(...)
#
#        bgColor: vector of any of the three kinds of R color specifications, 
#                 i.e., either a color name (as listed by colors()), or a 
#                 hexadecimal string of the form "#rrggbb" or "#rrggbbaa" 
#                 (see rgb), or a positive integer i meaning palette()[i].
#
#    Returned value:    None
#
#    Example:   drawBioNetNodeBackground(plotLocations)
#
#    Last revised on July 15, 2014
#

drawBioNetNodeBackground <- function(trackLocations, 
                        bgColor=gray(0.9, alpha=0.5)) {

    rgb.val <- as.vector(col2rgb(bgColor))/255
    bgColorColor <- rgb(red=rgb.val[1], green=rgb.val[2], 
                            blue=rgb.val[3], alpha=0.5)

    end <- nrow(trackLocations$inPositions)

    polygonX <- c(trackLocations$outPositions[1:end,1], 
                    trackLocations$inPositions[end:1,1])
    polygonY <- c(trackLocations$outPositions[1:end,2], 
                    trackLocations$inPositions[end:1,2])

    polygonX <- polygonX + trackLocations$nodeCenter[1]
    polygonY <- polygonY + trackLocations$nodeCenter[2]

    polygon(polygonX, polygonY, col=bgColor, border=NA)
}






#    __________________________________________________________________________
#    **************************************************************************
#
#    Draw quadratic Bezier curve between two samples inside a node. This is  
#    always in most inside of the node. 
#
#    Arguments:
#
#        nodeIndex: non negative integer, the node on which link line is drawn
#        fromSample: non negative integer, the first sample to be linked
#        toSample:   non negative integer, the second sample to be linked
#        outer:      non negative numeric, the start and end of link line  
#                    relative to node center
#        plotColors: character vector of R color specification
#
#    Return value:  None
#
#    Example:   linkBioNetSamples(nodeCenter,  fromSample=5,  toSample=90,
#                        outer, plotColors=red(1.0, alpha=0.5))
#
#    Last revised on July 15, 2014
#

linkBioNetSamples <- function(nodeIndex, fromSample, toSample, outer, 
                                    plotColors) {

    if(is.numeric(nodeIndex)  ==  FALSE || is.numeric(fromSample)  ==  FALSE ||
            fromSample<0 || is.numeric(toSample) == FALSE || toSample<0 || 
            is.numeric(outer) == FALSE || outer<0 || nodeIndex<0)
                stop("First four arguments must be non-negative integer.\n") 

    totalSample <- getBioNetPlotTotalSample()
    if(fromSample>totalSample || toSample>totalSample || 
                fromSample == toSample)
            stop("Incorrect sample index.\n")

    bioNetGraph <- getBioNetGraph()
    nodeRadius <- getBioNetNodeRadius()
    nodeCenter <- bioNetGraph$layout[nodeIndex, ]

    basePosition <- getBioNetBasePositions()
    sampleLocation <- getBioNetSamplePlotPosition(totalSample)
    pointLocations <- basePosition[sampleLocation[,2], 1:2]*outer*nodeRadius

    from <- as.numeric(pointLocations[fromSample, ])
    to <- as.numeric(pointLocations[toSample, ])
    linkLine <- getBezierCurve(from, to, 1000)

    lines(linkLine$posX+nodeCenter[1], linkLine$posY+nodeCenter[2], 
                col=plotColors)
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    Calculate x and y coordinated for a quadratic Bezier curve between two  
#    points with the equation:  
#
#        B(t) = (1-t) ((1-t)P0 + tP1) + t((1-t)P1 + tP2)
#
#    where P0 is the start point, P2 is the end point, and P1 is the control  
#    point. P1 will be adjusted based on the distance of two points.
#
#    Arguments:
#
#        lineStart: The point where Bezier line starts
#        lineEnd:   The point where Bezier line ends
#        totalPoints:   total number of points that form a circle line
#
#    Return value: a list contianing x and y coordinates for a quadratic 
#                  Bezier curve
#
#    Example:   internal use only
#
#    Last revised on July 15, 2014
#
#

getBezierCurve <- function(lineStart, lineEnd, totalPoints) {

    P0 <- as.numeric(lineStart)
    P2 <- as.numeric(lineEnd)

    t <- seq(0, 1, 1/totalPoints)
    linkX <- (1-t)^2*P0[1] + t^2*P2[1]
    linkY <- (1-t)^2*P0[2] + t^2*P2[2]

    return (list(posX=linkX, posY=linkY))
}




#    __________________________________________________________________________
#    **************************************************************************
#
#    Draw customized arrow between two nodes. 
#
#    Arguments:
#
#        fromNode:  non negative integer, the start node to be linked
#        toNode:    non negative integer, the end node to be linked
#        lineColor: character vector, color of the link line
#        arrowSize: non-negative numeric, scaling factor for arrow size, 
#                   default 1
#
#    Return value:  None
#
#    Example: linkBioNetNodes(from=5, to=10, lineColor="red", arrowSize=0.5)
#
#    Last revised on July 23, 2014
#

linkBioNetNodes <- function(fromNode, toNode, lineColor="black", arrowSize=1) {

    if(is.numeric(fromNode) == FALSE || is.numeric(toNode) == FALSE || 
            fromNode<0 || toNode<0 || fromNode == toNode)
        stop("Incorrect node number(s) defined.\n") 

    bioNetGraph <- getBioNetGraph()
    nodeLayout <- bioNetGraph$layout

    if(fromNode>nrow(nodeLayout) || toNode>nrow(nodeLayout))  
        stop("Incorrect node number(s) defined.\n") 

    fromCenter <- as.numeric(nodeLayout[fromNode, ])
    toCenter   <- as.numeric(nodeLayout[toNode, ])
    center2center <- sqrt((fromCenter[1]-toCenter[1])^2 + 
            (fromCenter[2]-toCenter[2])^2)

    basePositions <- getBioNetBasePositions()
    totalPoints <- ceiling((center2center/2/pi)*nrow(basePositions))

    lineX <- seq(fromCenter[1], toCenter[1], length=totalPoints)
    lineY <- seq(fromCenter[2], toCenter[2], length=totalPoints)

    nodeRadius <- getBioNetNodeRadius()
    plotAreaWid <- getBioNetPlotAreaWidth()
    plotRadius  <- nodeRadius*(1.5 + plotAreaWid)

    lineAdjust <- ceiling(plotRadius/center2center*totalPoints)
    lineX <- lineX[(lineAdjust+1):(totalPoints-lineAdjust)]
    lineY <- lineY[(lineAdjust+1):(totalPoints-lineAdjust)]

    last <- length(lineX)
    lineLength <- sqrt((lineX[1]-lineX[last])^2 + (lineY[2]-lineY[last])^2)

    if(lineLength<(arrowSize*nodeRadius)) 
        stop("No enough space for an arraow. Please reduce arrow size.\n") 

    nodeLinkLine <- getBioNetNodeLinkLine(lineX, lineY, arrowSize, lineLength)

    rgb.val <- as.vector(col2rgb(lineColor))/255
    lineColor <- rgb(red=rgb.val[1], green=rgb.val[2], blue=rgb.val[3])

    polygon(nodeLinkLine$arrowX, nodeLinkLine$arrowY, 
                    col=lineColor, border="black")
}





#    __________________________________________________________________________
#    **************************************************************************
#
#    Generate x and y coordinates for arrow head and tail with defined length. 
#    By default, the arrow is in inside of a circle (radius 1) without tail  
#    and it points to radian 0. The tail, if any, will be added to the left.
#
#    Arguments:
#
#        lineX:      numeric vector, x coordinates of the link line
#        lineY:      numeric vector, y coordinates of the link line
#        arrowSize:  non-negative numeric, scaling factor for arrow size, 
#                    default 1
#        lineLength: non negative integer, the length of link line
#
#    Return value:   x and y coordinates for link line with arrow
#
#    Example: line.corr <- getBioNetNodeLinkLine(lineX, lineY, arrowSize=1, 
#                               lineLength)
#    Last revised on July 23, 2014
#

getBioNetNodeLinkLine <- function(lineX, lineY, arrowSize=1, lineLength) {

    nodeRadius <- getBioNetNodeRadius()

    polygonX <- c(1,-0.7,-0.4,-1,  -1,  -0.4,-0.7,1)*nodeRadius
    polygonY <- c(0, 0.7, 0.2, 0.2,-0.2,-0.2,-0.7,0)*nodeRadius

    polygonX <- polygonX*arrowSize
    polygonY <- polygonY*arrowSize

    headLength <- max(polygonX)-min(polygonX)
    if(lineLength>headLength) {
        tailLength <- lineLength - headLength

        polygonX[4] <- polygonX[4] - tailLength
        polygonX[5] <- polygonX[5] - tailLength
        polygonX <- polygonX + (tailLength/2)

    } else {
        scaleFactor <- lineLength/headLength
        polygonX <- polygonX*scaleFactor
        polygonY <- polygonY*scaleFactor
    }

    mid <- ceiling(length(lineX)/2)
    lineCenter <- c(lineX[mid], lineY[mid])

    last <- length(lineX)
    angle <- atan2(lineY[last]-lineCenter[2], lineX[last]-lineCenter[1])

    newX <- polygonX
    newY <- polygonY

    for(a.point in seq_len(length(newX))) {
        newX[a.point] <- polygonX[a.point]*cos(angle) -  
                            polygonY[a.point]*sin(angle)
        newY[a.point] <- polygonX[a.point]*sin(angle) +  
                            polygonY[a.point]*cos(angle)
    }

    newX <- newX + lineCenter[1]
    newY <- newY + lineCenter[2]

    return (list(arrowX=newX, arrowY=newY))
}


#    Last revised on May 18, 2015
#    ________________________________________________________________________









