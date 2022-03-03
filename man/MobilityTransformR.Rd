\name{MobilityTransformR}
\alias{MobilityTransformR}
\docType{package}
\title{Mobility Transformation of CE-MS Data}
\description{
MobilityTransformR collects a tool set for effective mobility scale transformation of CE-MS data.
}
\details{
MobilityTransformR collects a tool set for effective mobility scale transformation of CE-MS/MS data in order to increase reproducibility. It provides functionality to determine the migration times from mobility markers that have been added to the analysis and performs the transformation based on these markers. MobilityTransformR supports the conversion of numeric vectors, Spectra-objects, and MSnOnDiskExp.
}
\author{
Author: \packageAuthor{MobilityTransformR}
Maintainer: \packageMaintainer{MobilityTransformR}
}
\references{
Breitling, R. et al. Ab initio prediction of metabolic networks using 
Fourier transform mass spectrometry data. 2006. Metabolomics 2: 155--164. 
10.1007/s11306-006-0029-z
}
\keyword{mass spectrometry, metabolomics}

\examples{
marker <- data.frame(markerID = c("marker1", "marker2"),
                      rtime = c(20,80),
                      mobility = c(0, 2000))
                      
## transform vector of migration times                       
rtime <- c(10,20,30,40,50,60,70,80,90,100)
mobilityTransform(x = rtime, marker = marker)
                       
}
