## ----style, echo = FALSE, results = 'asis', message=FALSE---------------------
BiocStyle::markdown()

## ---- echo = FALSE, message = FALSE-------------------------------------------
library(BiocStyle)

## ---- message=FALSE, warning=FALSE, message=FALSE-----------------------------
devtools::install_github("LiesaSalzer/MobilityTransformR")

## ----libraries, message=FALSE, warning=FALSE----------------------------------
# load required libraries
library(MobilityTransformR)
library(xcms)
library(Spectra)

## ----data, message=FALSE------------------------------------------------------
fl <- list.files(system.file("extdata/", package = "MobilityTransformR"), 
                 pattern = ".mzML", full.names = T)

# Load mzXML data with MSnBase
raw_data <- readMSData(files = fl, mode = "onDisk")

## ----EOF marker, message=FALSE, warning=FALSE---------------------------------
# mz tolerance depends on MS mass accuracy 
tolerance <- 0.005

# [M+H]+ of paracetamol: mz = 152.071154
mz_paracetamol <- c(152.071154 - tolerance, 152.071154 + tolerance)
mt_paracetamol <- c(600, 1500)

marker_EIE <-  raw_data %>% filterMz(mz = mz_paracetamol) %>% 
    filterRt(rt = mt_paracetamol)

plot(chromatogram(marker_EIE), 
     main = "Paracetamol EIE", 
     xlab = "migration time (sec)")

# adjust mz and MT windows if necessary

## ----EOF marker getMT, message=FALSE, warning=FALSE---------------------------
# get the MT of paracetamol
paracetamol <- getMtime(raw_data,
                        mz = mz_paracetamol, 
                        mt = mt_paracetamol)
paracetamol


## ----charged marker, message=FALSE, warning=FALSE-----------------------------
mz_procaine <- c(237.160303 - tolerance, 237.160303 + tolerance)
mt_procaine <- c(300, 800)

marker_EIE <-  raw_data %>% filterMz(mz = mz_procaine) %>% 
    filterRt(rt = mt_procaine)

plot(chromatogram(marker_EIE), 
     main = "Procaine EIE",
     xlab = "migration time (sec)")

# get the MT of procaine using adjusted parameters

procaine <- getMtime(raw_data,
                        mz = mz_procaine, 
                        mt = mt_procaine)

procaine

## ---- message=FALSE-----------------------------------------------------------
# Create a data.frame that stores marker information 
marker <- paracetamol
marker$markerID = "Paracetamol"
marker$mobility = 0

## ----mobility, message=FALSE--------------------------------------------------
procaineMobility <- mobilityTransform(x = procaine[1,1], marker = marker[1,], 
                                      tR = 3/60, U = +30, L = 800)
procaineMobility


## ---- message=FALSE-----------------------------------------------------------
procaine$markerID = "Procaine"
procaine$mobility = procaineMobility
marker <- rbind(marker, procaine)

## ---- message=FALSE-----------------------------------------------------------
# Conversion of mt in OnDiskMSnExp objects
mobility_data <- mobilityTransform(x = raw_data, marker = marker)

# #OnDiskMSnExp can by exported by writeMSData, Note that it is important to set
# #copy = FALSE, (otherwise spectrum ordering will be wrong)
# fl_mobility_data <- tempfile()
# writeMSData(filterFile(mobility_data, 1), file = fl_mobility_data, copy = FALSE)

## ---- message=FALSE-----------------------------------------------------------
#load the test data as spectra object
spectra_data <- Spectra(fl[1], backend = MsBackendMzR())

spectra_mobility <- mobilityTransform(spectra_data, 
                                      marker[marker$fileIdx == 1,])


# #Transformed data can then be exported again as .mzML file to use in xcms or 
# #other software
# fl_mobility <- tempfile()
# export(spectra_mobility, MsBackendMzR(), file = fl_mobility)


## ----mobility transformed data, message=FALSE---------------------------------
# Example: Extract ion electropherogram (EIE) from lysine
mz_lysine <- c(147.112806
 - tolerance, 147.112806
 + tolerance)
mobilityRestriction <- c(-200, 2500)


# Extract ion electropherogram of compound
lysine_EIE <-  mobility_data %>% 
    filterMz(mz = mz_lysine) %>% 
    filterRt(rt = mobilityRestriction)

plot(chromatogram(lysine_EIE), 
     main = expression(paste("Lysine EIE - µ"[eff]," scale")), 
     xlab = expression(paste("µ"[eff],"  (", frac("mm"^2,"kV min"),")"))
     )

# compare with extracted ion electropherogram of migration time scale
lysine_mt_EIE <-  raw_data %>% 
    filterMz(mz = mz_lysine) 

plot(chromatogram(lysine_mt_EIE),
     main = "Lysine EIE - MT scale", 
     xlab = "MT (sec)")


## ----si-----------------------------------------------------------------------
sessionInfo()

