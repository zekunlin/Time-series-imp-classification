getwd()
data_repo <- "D:/Zekun/Landsat_ARD/Landsat_ARD/wake_ard_period2/"
data_repo <- "D:/Zekun/Landsat_ARD/Landsat_ARD/wake_ard_period1/test/"
list.files(data_repo, pattern=".part")
library(raster)
B1 <- list.files(data_repo, pattern = "SRB1", full.names = T)
B2 <- list.files(data_repo, pattern = "SRB2", full.names = T)
B3 <- list.files(data_repo, pattern = "SRB3", full.names = T)
B4 <- list.files(data_repo, pattern = "SRB4", full.names = T)
B5 <- list.files(data_repo, pattern = "SRB5", full.names = T)
B7 <- list.files(data_repo, pattern = "SRB7", full.names = T)
BQA <- list.files(data_repo, pattern = "PIXELQA_doy", full.names = T)



for (i in B2) {
    bandx <- raster(B1[i])
}
b1 <- raster(B1[1])
plot(b1)

