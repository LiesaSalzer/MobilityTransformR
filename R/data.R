#' @name data
#'
#' @aliases data
#'
#' @title CE-MS test data
#'
#' @description
#' The test data of the `MobilityTransformR` package consist of two files,
#' i.e. `"CEMS_10ppm.mzML"` and `"CEMS_25ppm.mzML"`.
#' The data contains CE-MS runs of a standard mixture that contains e.g.
#' Lysine (at 10 ppm and 25 ppm, respectively) and the neutral EOF marker
#' Paracetamol (50 ppm).
#' The data was acquired on a 7100 capillary electrophoresis system from
#' Agilent Technologies, coupled to an Agilent 6560 IM-QToF-MS.
#' CE Separation was performed using a 80 cm fused silica capillary with an
#' internal diameter of 50 µm and external diameter of 365 µm.
#' The Background Electrolyte was 10 % acetic acid and separation was
#' performed at +30 kV and a constant pressure of 50 mbar.
#' MS detection was performed in positive ionization mode.
#'
#' The raw data were then converted to an open-source ".mzML" format
#' (Proteowizzard) and load into R via the `MSnBase::readMSData()` function.
#' In order to reduce data size, the test data was subsequently cutted in
#' migration time and mz range using `filterRt(rt = c(400, 900))`
#' and `filterMz(mz = c(147.1, 152.0))` from `MSnBase`
#'
#' @author Liesa Salzer
#'
