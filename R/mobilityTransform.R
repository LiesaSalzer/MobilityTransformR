#' @name mobilityTransform
#'
#' @aliases mobilityTransform
#'
#' @title Effective mobility scale transformation
#'
#' @description
#' The function `mobilityTransform` performs effective mobility scale 
#' transformation of CE data 
#'
#' @param tM
#' `numeric` migration time `tM` in sec which will be transformed to an 
#' effective mobility scale
#'
#' @param type
#' `character` defines if one or two markers are beeing used for effective 
#' mobility scale transformation. 
#' If `type` is `singleMarker` either `tEOF` or `tA` (in combination with 
#' `mobilityA`) must be provided along with `U`, `l`, and `L`
#' If `type` is `twoMarker`, `tEOF`, `tA`, and `mobilityA` must be provided.
#' 
#'  @param tEOF
#' `numeric` migration time of the (neutral charged) EOF marker in sec
#' 
#'  @param tA
#' `numeric` migration time of the charged EOF marker in sec
#' 
#'  @param tR
#' `numeric` field ramping delay in sec, default is 0
#' 
#'  @param mobilityA
#' `numeric` effective mobility of the charged marker `tA` in mm^2 / (kV * min)
#' 
#'  @param U
#' `numeric` applied electrical voltage U in kV. Note that in reversed CE 
#' polarity, also a negative voltage is applied.
#' 
#'  @param L
#' `numeric` total length of the capillary in mm
#' 
#'  @param l
#' `numeric` effective length of the capillary in mm. In CE-MS l = L and in 
#' CE-UV l < L. The default is l = L.
#'
#' @details
#' `mobilityTransform` performs effective mobility scale transformation of 
#' time-scaled CE data `tM`. 
#' For the transformation, either one mobility marker (`type` = "singleMarker") 
#' or both (`type` = "towMarker") were used. 
#' If a single marker is used, either a neutral EOF marker `tEOF`, or charged 
#' marker `tA` with its corresponding mobility `mobilityA` must be provided, 
#' along with the applied voltage `U`, and the total and effective capillary 
#' length `L` and `l`.
#' If two markers are used, both a neutral EOF marker `tEOF`, or charged 
#' marker `tA` with its corresponding mobility `mobilityA` must be provided.
#' Additionally, field ramping delays can be included by `tR`, which will 
#' result in more precise effective mobility values.  
#' 
#' @return
#' `numeric` that represents effective mobility in mm^2 / (kV * min)
#'
#' @author Liesa Salzer, \email{liesa.salzer@@helmholtz-muenchen.de}
#'
#' @examples
#'  ######  example to be added!
#'
#'@export
mobilityTransform <- function(tM, type = c("singleMarker", "twoMarker"),
                              tEOF = NA, tA = NA, tR = 0, mobilityA = NA,
                              U = NA, L = NA, l = L,  ...) {
  
  if (!any(c("singleMarker", "twoMarker") %in% type)) {
    stop("'type' must contain either 'singleMarker' or 'twoMarker'")}
  
  if (missing(tM)) {stop("'tM' is missing")}
  
  tM <- tM/60
  tR <- tR/60
  
  if ("singleMarker" %in% type) {
    
    if (is.na(tEOF) & is.na(tA)) {
      stop("Either 'tEOF' or 'tA' must be provided")}
    if (missing(U)) {stop("'U' is missing")}
    if (missing(L)) {stop("'L' is missing")}
    
    E <- U / L
    
    
    
    if (!is.na(tEOF)) {
      tEOF <- tEOF/60
      mobility <- l / E * ((1 / (tM - (tR / 2))) - (1 / (tEOF - (tR / 2))))
    }
    
    if (!is.na(tA)) {
      tA <- tA/60
      if (is.na(mobilityA)) {stop("'mobilityA' is missing")}
      
      mobility <- mobilityA + l / E * ((1 / (tM - (tR / 2))) - (1 / (tA - (tR / 2))))
    }
    
  } else if ("twoMarker" %in% type) {
    
    if (is.na(tEOF) | is.na(tA)) {
      stop("Both 'tEOF' and 'tA' must be provided")}
    
    if (is.na(mobilityA)) {stop("'mobilityA' is missing")}
    tEOF <- tEOF/60
    tA <- tA/60
    mobility <- ((tM - tEOF) /
                   (tA - tEOF)) * ((tA - (tR / 2)) /
                                     (tM - (tR / 2))) * mobilityA
  }
  
  
  return(mobility)
}

