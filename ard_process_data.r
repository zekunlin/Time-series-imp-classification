getwd()
#data_repo <- "D:/Zekun/Landsat_ARD/Landsat_ARD/wake_ard_period2/"
data_repo <- "D:/Zekun/Landsat_ARD/Landsat_ARD/wake_ard_period2/"

library(raster)
library(stringr)
library(matrixStats)

## now do the image composite by seasons based on DOY

img_date <- str_extract(b1_11, '[0-9][0-9][0-9][0-9][0-9][0-9][0-9]') # get doy of each image
yr <- 2011
origins<- paste(yr, "-01-01", sep = "")
substr(img_date, 5, 7)
img_mth_days <- as.Date(as.numeric(substr(img_date, 5,7)), origin = origins)
mth_numbers <- as.numeric(format(img_mth_days, "%m"))

s1_ind <- which(mth_numbers <= 3)
s2_ind <- which(mth_numbers <= 6 & mth_numbers > 3)
s3_ind <- which(mth_numbers <= 8 & mth_numbers > 6)
s4_ind <- which(mth_numbers <= 12 & mth_numbers >8)
seasons_ind <- list(s1 = s1_ind, s2 = s2_ind, s3 = s3_ind, s4 = s4_ind)
##############################################################################################

## now we are getting files for all sr data
b1 <- list.files(data_repo, pattern = "SRB1", full.names = T)
b2 <- list.files(data_repo, pattern = "SRB2", full.names = T)
b3 <- list.files(data_repo, pattern = "SRB3", full.names = T)
b4 <- list.files(data_repo, pattern = "SRB4", full.names = T)
b5 <- list.files(data_repo, pattern = "SRB5", full.names = T)
b7 <- list.files(data_repo, pattern = "SRB7", full.names = T)
bqa <- list.files(data_repo, pattern = "PIXELQA_doy", full.names = T)

b1_11 <- list.files(data_repo, pattern = "^.*SRB1.doy2011.*", full.names = T)
b2_11 <- list.files(data_repo, pattern = "^.*SRB2.doy2011.*", full.names = T)
b3_11 <- list.files(data_repo, pattern = "^.*SRB3.doy2011.*", full.names = T)
b4_11 <- list.files(data_repo, pattern = "^.*SRB4.doy2011.*", full.names = T)
b5_11 <- list.files(data_repo, pattern = "^.*SRB5.doy2011.*", full.names = T)
b7_11 <- list.files(data_repo, pattern = "^.*SRB7.doy2011.*", full.names = T)
bqa_11 <- list.files(data_repo, pattern = "^.*PIXELQA.doy2011.*", full.names = T)

# stack and convert to matrix
b1_stack <- stack(b1_11)
b1_matrix <- as.matrix(b1_stack)
b2_stack <- stack(b2_11)
b2_matrix <- as.matrix(b2_stack)
b3_stack <- stack(b3_11)
b3_matrix <- as.matrix(b3_stack)
b4_stack <- stack(b4_11)
b4_matrix <- as.matrix(b4_stack)
b5_stack <- stack(b5_11)
b5_matrix <- as.matrix(b5_stack)
b7_stack <- stack(b7_11)
b7_matrix <- as.matrix(b7_stack)

b1_matrix_copy <- b1_matrix
bqa_stack <- stack(bqa_11)
qa_matrix <- as.matrix(bqa_stack)

bad_pixel_flag <- c(1, 68, 72, 80, 96, 132, 136, 144, 160, 224)
n_layer <- ncol(b1_matrix)

for(y in 1:n_layer){
    # we want the clear pixels only, so remove other pixels keep flag the 66 and 130 and do this for all bands
    b1_matrix[, y] <- replace(b1_matrix[,y], qa_matrix[,y] %in% bad_pixel_flag, NA)  
    b2_matrix[, y] <- replace(b2_matrix[,y], qa_matrix[,y] %in% bad_pixel_flag, NA)
    b3_matrix[, y] <- replace(b3_matrix[,y], qa_matrix[,y] %in% bad_pixel_flag, NA)
    b4_matrix[, y] <- replace(b4_matrix[,y], qa_matrix[,y] %in% bad_pixel_flag, NA)
    b5_matrix[, y] <- replace(b5_matrix[,y], qa_matrix[,y] %in% bad_pixel_flag, NA)
    b7_matrix[, y] <- replace(b7_matrix[,y], qa_matrix[,y] %in% bad_pixel_flag, NA)
}



sr_allBands_matrix <- list(band1 = b1_matrix, band2 = b2_matrix, band3 = b3_matrix, 
                            band4 = b4_matrix, band5 = b5_matrix, band7 = b7_matrix)
ndvi <- function(band3, band4){
    ndvi_val <- (band4 - band3)/(band4 + band3)
    return(ndvi_val)
}

evi <- function(band1, band3, band4){
    evi_val <- (2.5 * ((band4 - band3)/(band4 + 6*band3 - 7.5*band1 +1)))
    return(evi_val)
}

ndvi_matrix <- (sr_allBands_matrix$band4 - sr_allBands_matrix$band3)/(sr_allBands_matrix$band4 + sr_allBands_matrix$band3)
ndvi_by_season <- list(
    s1 = rowMedians(ndvi_matrix[,s1_ind]),
    s2 = rowMedians(ndvi_matrix[,s2_ind]),
    s3 = rowMedians(ndvi_matrix[,s3_ind]),
    s4 = rowMedians(ndvi_matrix[,s4_ind])
)

s1.bandx <- lapply(sr_allBands_matrix, "[", ,seasons_ind[["s1"]])
s1.med <- lapply(s1.bandx, rowMedians, na.rm = TRUE)

s2.bandx <- lapply(sr_allBands_matrix, "[", ,seasons_ind[["s2"]])
s2.med <- lapply(s1.bandx, rowMedians, na.rm = TRUE)

s3.bandx <- lapply(sr_allBands_matrix, "[", ,seasons_ind[["s3"]])
s3.med <- lapply(s1.bandx, rowMedians, na.rm = TRUE)

s4.bandx <- lapply(sr_allBands_matrix, "[", ,seasons_ind[["s4"]])
s4.med <- lapply(s1.bandx, rowMedians, na.rm = TRUE)



lapply(sr_allBands_matrix, FUN = matrixStats::rowMedians, na.rm = TRUE)
s1_b1 <- b1_matrix[, s1_ind]
s1_b1_med <- matrixStats::rowMedians(s1_b1, na.rm = TRUE, dim. = dim(s1_b1))

