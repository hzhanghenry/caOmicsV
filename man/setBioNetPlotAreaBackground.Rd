\name{setBioNetPlotAreaBackground}
\alias{setBioNetPlotAreaBackground}
\title{
    Set up Plot Area Background for Nodes of caOmicsV bioNetCircos Layout
}
\description{
    Change the plot area background of igraph nodes. Use white color to erase 
    background and grey to show the plot area boundary.
}
\usage{
    setBioNetPlotAreaBackground(bgColor=grey(0.75, alpha=0.5))
}
\arguments{
    \item{bgColor}{
    character vector of a colour name or a R colour specification
    }
}
\value{
    None
}
\author{
    Henry Zhang
}
\examples{
    setBioNetPlotAreaBackground(grey(0.75, alpha=0.5))
}
\keyword{ methods }
