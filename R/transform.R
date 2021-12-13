#' @name mobilityTransform
#'
#' @aliases mobilityTransform
#'
#' @title Effective mobility scale transformation of CE-MS data
#'
#' @description
#' `mobilityTransform` performs effective mobility scale transformation of CE(-MS) 
#' data, which is used to overcome variations of the migration times, caused by 
#' differences in the Electroosmotic Flow (EOF) between different runs. 
#' In order to monitor the EOF and perform the transformation, neutral or 
#' charged EOF markers are spiked into the sample before analysis. The 
#' information of the EOF markers (migration time and effective mobility) will 
#' be then used to perform the  effective mobility transformation of the 
#' migration time scale. 
#' 
#' For the transformation, either one mobility marker or both can be used. 
#' If a single marker is used, either a neutral EOF marker, or charged 
#' marker with its corresponding mobilities (0 for the neutral marker) must be
#' provided, along with the applied voltage `U`, and the total capillary length 
#' `L`.
#' If two markers are used, both a neutral EOF marker and a charged marker 
#' including their corresponding mobility must be provided.
#' Additionally, field ramping delays can be included by `tR`, which will 
#' result in more precise effective mobility values.  
#' 
#' Currently, `mobilityTransform` supports `numeric` vectors of migration times 
#' as input, `Spectra`-objects or `MSnOnDiskExp`-objects. 
#' `mobilityTransform` is a method that used different functions to convert 
#' CE-MS data, depending on the input class. Following functions will be applied 
#' depending on the input class:
#' 
#'  - `.transformNumeric`: performs effective mobility scale transformation of a 
#'    `numeric` migration time vector as input. This can be used to transform a 
#'    row of migration times or a single value. This function will return a 
#'    `numeric`
#' 
#' 
#'  - `.transformSpectra` performs effective mobility scale transformation of 
#'    the migration time scale within a `Spectra` object. This function will 
#'    return a `Spectra` object with effective mobility scale
#' 
#' 
#'  - `.transformOnDiskMSnExp` performs effective mobility scale transformation 
#'    of the migration time scale within an `OnDiskMSnExp` object. 
#'    Since `OnDiskMSnExp` can store multiple files, it is also possible to 
#'    perform the transformation of multiple files. Hence, 
#'    `.transformOnDiskMSnExp` requires the `marker` `data.frame` to have an 
#'    additional column "fileIdx" that stores the file Index of the migration 
#'    time of all markers. 
#' 
#' @param x
#' `numeric` migration time vector, `Spectra` -object, or `MSnOnDiskExp`-object 
#' that serves as input file to perform the effective mobility transformation.
#' The respective migration time scale should be in seconds(!).
#'
#' @param marker
#' `data.frame` containing minimum two columns, where one holds the determined 
#' migration time in minutes (here referred to as "rtime") of the EOF marker in 
#' the same run in which the migration time is going to be transformed and the 
#' other column the respective mobility ("mobility") of the EOF markers. Each 
#' row hold the values for one EOF marker. 
#' If `OnDiskMSnExp` is used in `x`, a third column "fileIdx" is required, that
#' stores the file Index. 
#' One or two entries are required per file for the transformation and depending 
#' on the number of entries the transformation will be performed either on one or 
#' both markers.
#' 
#'@param ...
#' Additional parameters passed to `convertMtime`, as for example `L`, `U`, 
#' and `tR`.
#'  
#' 
#' @return
#' The same class as the input class will be returned, i.e. if a `numeric` is 
#' used as input a `numeric` that represents effective mobility will be returned.
#' If a `Spectra`-Object is the input, also a `Spectra`-Object with transformed 
#' mobility scale will be returned. The same applies for `MSnOnDiskExp`-objects. 
#' The respective unit for the effective mobility is mm^2 / (kV * min) 
#'
#' @author Liesa Salzer
#'
#' @importFrom MetaboCoreUtils convertMtime
#'
#' @examples
#'  rtime <- c(10,20,30,40,50,60,70,80,90,100)
#'  marker <- data.frame(markerID = c("marker1", "marker2"),
#'                       rtime = c(20,80),
#'                       mobility = c(0, 2000))
#'  mobilityTransform(x = rtime, marker = marker)
#'
#'@export

mobilityTransform <- function(x, marker, 
                              tR = 0, 
                              U = numeric(), 
                              L = numeric()) {
  
  ## sanity checks
  if (!class(x) %in% c("numeric","Spectra","OnDiskMSnExp")) 
    stop("'x' needs to be of class 'numeric', 'Spectra' or 'OnDistMSnExp' but 
          not class '", class(x),"'")
  if (missing(marker)) {
    stop("Missing data.frame 'marker' with marker information")}
  if (!all(c("rtime","mobility") %in% colnames(marker))) {
    stop("Missing column 'rtime', 'mobility' or both")}
  
  
  if (is(x, "numeric")) {
    FUN = .transformNumeric
    FUN <- match.fun(FUN)
    do.call(FUN, list(x = x, marker = marker, tR = tR, U = U, L = L))
  }
  else if (is(x, "Spectra")) {
    FUN = .transformSpectra
    FUN <- match.fun(FUN)
    do.call(FUN, list(x = x, marker = marker, tR = tR, U = U, L = L))
  }
  
  else if (is(x, "OnDiskMSnExp")) {
    FUN = .transformOnDiskMSnExp
    FUN <- match.fun(FUN)
    do.call(FUN, list(x = x, marker = marker, tR = tR, U = U, L = L))
  }
  
  
}

#' @name .transformNumeric
#' 
#' @title Transformation of numeric
#'
#' @param x
#' `numeric` migration time vector in seconds.
#' 
#' @return
#' `numeric` vector that represents effective mobility in mm^2 / (kV * min) 
#' 
#' @example 
#' rtime <- c(10,20,30,40,50,60,70,80,90,100)
#' marker <- data.frame(markerID = c("marker1", "marker2"),
#'                      rtime = c(20,80),
#'                      mobility = c(0, 2000))
#' transformNumeric(x = rtime, marker = marker)
#' transformNumeric(x = rtime, marker = marker[-1,], U = 30, L = 90)
#' 
.transformNumeric <- function(x, marker, tR = tR, U = U, L = L) {

  convertMtime(x/60, 
               rtime = marker$rtime/60, 
               mobility = marker$mobility, tR = tR, U = U, L = L)
  
}

#' @name .transformSpectra
#' 
#' @title Transformation of Spectra
#'
#' @param x
#' `Spectra`-object that stores the migration times in seconds.
#' 
#' @return
#' `Spectra`-Object that stores the effective mobility in mm^2 / (kV * min). 
#' 
#' @example 
#' fl <- system.file("extdata/CEMS_metabolites_10ppm_pos_centroidedData.mzML", 
#' package = "MobilityTransformationR")
#' spectra_data <- Spectra(fl, backend = MsBackendMzR())
#' marker <- data.frame(markerID = c("marker1", "marker2"),
#'                      rtime = c(20,80),
#'                      mobility = c(0, 2000))
#' .transformSpectra(x = spectra_data, marker = marker)

.transformSpectra <- function(x, marker, tR = tR, U = U, L = L) {
  
  xTransf <- x
  xTransf$rtime <- convertMtime(xTransf$rtime/60, 
                            rtime = marker$rtime/60, 
                            mobility = marker$mobility, tR = tR, U = U, L = L)
  
  ## Data needs to be ordered by the migration time and spectrum IDs needs to 
  ## be removed to prevent errors in xcms 
  xTransf <- xTransf[order(rtime(xTransf))]
  xTransf$spectrumId <- NA
  
  return(xTransf)
  
}

#' @name .transformOnDiskMSnExp
#' 
#' @title Transformation of OnDiskMSnExp
#' 
#' @param x
#' `OnDiskMSnExp`-object that stores the migration times in seconds.
#'  
#' @return
#' `OnDiskMSnExp`-Object that stores the effective mobility in mm^2 / (kV * min). 
#' 
#' 
#' @examples 
#' fl <- system.file("extdata/CEMS_metabolites_10ppm_pos_centroidedData.mzML", 
#' package = "MobilityTransformationR")
#' raw_data <- readMSData(files = fl,
#'                        mode = "onDisk")
#' marker <- data.frame(markerID = c("marker1", "marker2"),
#'                      rtime = c(20,80),
#'                      mobility = c(0, 2000),
#'                      fileIdx = c(1,1))
#'                      
#' .transformOnDiskMSnExp(x = raw_data, marker = marker)

.transformOnDiskMSnExp <- function(x, marker, tR = tR, U = U, L = L) {
  ## sanity checks
  if (!all(c("fileIdx") %in% colnames(marker))) {
    stop("Missing column 'fileIdx'")}
  
  xTransf <- x
  
  rt_file <- split(rtime(x), fromFile(x))
  
  for (i in names(rt_file)) {
    ## Filter fData(xTransf) based on file index and change rt into transformed scale
    fData(xTransf)[fData(xTransf)$fileIdx == i,]$retentionTime <- 
      convertMtime(x = rt_file[[i]]/60, 
                   rtime = marker[marker$fileIdx == i,]$rtime/60, 
                   mobility = marker[marker$fileIdx == i,]$mobility, tR = tR, U = U, L = L)
    
    ## Data needs to be ordered by the migration time and spectrum IDs needs to 
    ## be removed to prevent errors in xcms 
    fData(xTransf)[fData(xTransf)$fileIdx == i,]$retentionTime <- 
      order(fData(xTransf)[fData(xTransf)$fileIdx == i,]$retentionTime)
    
  }
  
  fData(xTransf)$spectrumId <- NA
  
  return(xTransf)
  
}
