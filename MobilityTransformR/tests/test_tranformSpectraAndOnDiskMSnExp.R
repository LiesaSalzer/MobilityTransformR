## Control tests 


test_data = readMSData(files = fl_mobility_data, mode = "onDisk")

test_data <- readMSData(files = fl_mobility, mode = "onDisk")

test_data = mobility_data



# mz tolerance depends on MS mass accuracy 
tolerance <- 0.005

# [M+H]+ of paracetamol: mz = 152.071154
mz_paracetamol <- c(152.071154 - tolerance, 152.071154 + tolerance)

EIE_paracetamol <-  test_data %>% filterMz(mz = mz_paracetamol) %>% filterRt(c(-200,1000))

plot(chromatogram(EIE_paracetamol), 
     main = "Paracetamol EIE", 
     sub = "Paracetamol expected @ 0 mobility",
     xlab = "migration time (sec)")




# align cwp 
cwp <- CentWaveParam(peakwidth = c(20, 150), noise = 1000, )
mfp <- MatchedFilterParam(binSize =1, snthresh = 3)
mq <- MassifquantParam(peakwidth = c(10,100), snthresh = 100, withWave = T)

xdata <- findChromPeaks(EIE_paracetamol, 
                        param = mfp)


# test if peaks were merged
plot(chromatogram(xdata))


c_peaks <- as.data.frame(chromPeaks(xdata))
c_peaks$CP <- row.names(c_peaks)
c_peaks$rtrange <- (c_peaks$rtmax-c_peaks$rtmin)


