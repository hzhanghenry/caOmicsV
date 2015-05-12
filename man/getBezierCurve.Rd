\name{getBezierCurve}
\alias{getBezierCurve}
\title{
    Calculate x and y Coordinates for A Quadratic Bezier Curve
}
\description{
    Calculate x and y coordinates for a quadratic Bezier curve between two 
    points with the equation:  
        B(t) = (1-t) ((1-t)P0 + tP1) + t((1-t)P1 + tP2)  
    where P0 is the start point, P2 is the end point, and P1 is the control 
    point. P1 will be adjusted based on the distance of two points.
}
\usage{
    getBezierCurve(lineStart, lineEnd, totalPoints)
}
\arguments{
    \item{lineStart}{
    numeric vector, the coordinate of a point where Bezier line starts
    }
    \item{lineEnd}{
    numeric vector, the coordinate of a point where Bezier line ends
    }
    \item{totalPoints}{
    non-negative numeric, total number of points that form a Bezier line
    }
}
\value{
    \item{posX}{x coordinates of points that form Bezier line}
    \item{posY}{y coordinates of points that form Bezier line}
}
\author{
    Henry Zhang
}
\examples{
    lineStart   <- c(0, 1)
    lineEnd     <- c(1, 0)
    totalPoints <- 2000
    the_line <- getBezierCurve(lineStart, lineEnd, totalPoints)
}
\keyword{ methods }

