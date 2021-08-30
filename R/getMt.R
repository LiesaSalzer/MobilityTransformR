#' @name getMt
#'
#' @aliases getMt
#'
#' @title get migration time of a known compound in CE-MS data
#'
#' @description
#' The function `getMt` searches the migration time (mt) of a known compound 
#' within a specified mz and mt range (`mz` and `mt`). 
#' Also, a minimum intensity value `minInt` of the
#' compound must be provided to filter from noise.
#'
#' @param x
#' `OnDiskMSnExp` object, containing CE-MS(/MS) data from a single file (e.g.
#' filtered by `filterFile()`). 
#' The migration time (mt) is provided in sec. The CE-MS data contain known 
#' compounds `mz`, of which the exact migration time is going to be determined. 
#'
#' @param mz
#' `numeric`, with upper and lower limit, (length(`mz`) = 2), 
#'  mz range of the compound of which the exact migration time 
#'  is going to be determined. The range should be as narrow as possible but 
#'  depends on the mass accuracy of the mass spectrometer that has been used 
#'  to acquire the data. 
#'
#' @param mt
#' `numeric`, with upper and lower limit, (length(`mt`) = 2),
#' mt range of where the compound is expected to migrate. 
#' Use a narrow `mt` range, to to avoid that other components with the same mz 
#' and different mt are determined and to speed up the calculations. 
#' 
#' @param minInt
#' `numeric`, minimum intensity of the compound of which the exact migration 
#' time is going to be determined. The `minInt` should be min at S/N = 3, to 
#' ensure that no noise is beeing picked
#' 
#'
#' @details
#' `getMt` uses `OnDiskMSnExp` object of CE-MS data to search for the 
#' migration time of selected compounds as for example EOF markers that where 
#' included to the analysis to later perform effective mobility scale 
#' transformation. 
#' The `OnDiskMSnExp` object is beeing filtered using the defined mz-range `mz`,
#' which depends on the mass accuracy of the mass spectrometer, and the 
#' -range `mt`, where the compound is expected to migrate. 
#' The migration time is determined by searching the mt of the maximum intense 
#' data point in the produced extracted ion electropherogram. Therefore,
#' `minInt` is used to filter from noise. 
#' 
#'
#' @return
#' single `numeric` that represents the migration time in sec 
#'
#' @author Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
#'
#' @examples
#'  ######  example to be added!
#'
#'@export
getMt <- function(x, mz, mt, minInt, ...) {
# sanity checks
if (!is(x, "OnDiskMSnExp")) {
  stop("'x' is not of class 'OnDiskMSnExp'!")
}
if (missing(mz) | missing(mt) | missing(minInt)) {
  stop(
    "Arguments 'mz', 'mt', and 'minInt'",
    " are required!"
  )
}
if (length(minInt) != 1) {
  stop("'minInt' needs to be of length 1")
}
if (!is.double(mz)) {
  mz <- as.double(mz)
}
if (!is.double(mt)) {
  mt <- as.double(mt)
}
if (!is.double(minInt)) {
  minInt <- as.double(minInt)
}


marker <- x %>%
  filterMz(mz = mz) %>%
  filterRt(rt = mt)

eic_marker <- chromatogram(marker, aggregationFun = "max")


df <- as.data.frame(intensity(eic_marker[, 1]))
df <- add_column(df, as.data.frame(rtime(eic_marker[, 1])))

maxInt <- df[which.max(df[, 1]), 1]
maxMt <- df[which.max(df[, 1]), 2]

if (maxInt < minInt) {
  stop("No peaks have been found. Align input parameters.")
}


return(maxMt)
}
