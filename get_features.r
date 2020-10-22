library(ncdf4)
library(stringr)
library(matrixStats)

bad_pixel_flag <- c(1, 68, 72, 80, 96, 132, 136, 144, 160, 224) # mask these pixels with flags from Landsat User Guide
## now do the image composite by seasons based on DOY


for (yr in 1984:2011) {
    ncf.in <- nc_open(paste("D:/Zekun/Imp_classification/ard_nc_files/LT05_Wake_ncdf4_", yr, ".nc", sep = ""))
    T <- ncvar_get(ncf.in, 'time')
    b1 <- ncvar_get(ncf.in, 'b1_sr')
    b2 <- ncvar_get(ncf.in, 'b2_sr')
    b3 <- ncvar_get(ncf.in, 'b3_sr')
    b4 <- ncvar_get(ncf.in, 'b4_sr')
    b5 <- ncvar_get(ncf.in, 'b5_sr')
    b7 <- ncvar_get(ncf.in, 'b7_sr')
    bqa <- ncvar_get(ncf.in, 'qa')
    lon.in <- ncvar_get(ncf.in, 'lon')
    lat.in <- ncvar_get(ncf.in, 'lat')
    crs.in <- ncatt_get(ncf.in, 0,'CRS')
    cell.in <- ncatt_get(ncf.in, 0,'cell_size')
    xmin.in <- ncatt_get(ncf.in, 0,'xmin')
    xmax.in <- ncatt_get(ncf.in, 0,'xmax')
    ymin.in <- ncatt_get(ncf.in, 0,'ymin')
    ymax.in <- ncatt_get(ncf.in, 0,'ymax')

    nc_close(ncf.in)
    
    nc.ncol <- dim(b1)[2]
    nc.nrow <- dim(b1)[1]
    img_doy <- str_extract(T, '[0-9][0-9][0-9][0-9][0-9][0-9][0-9]') # get doy of each image
    origins <- paste(yr, "-01-01", sep = "")
    img_date <- as.Date(as.numeric(substr(img_doy, 5,7)), origin = origins)
    mth <- as.numeric(format(img_date, "%m"))

    n_layer <- dim(b1)[3]
    for(y in 1:n_layer){
    # we want the clear pixels only, so remove other pixels keep flag the 66 and 130 and do this for all bands
        b1[,,y] <- replace(b1[,,y], bqa[,,y] %in% bad_pixel_flag, NA)  
        b2[,,y] <- replace(b2[,,y], bqa[,,y] %in% bad_pixel_flag, NA)
        b3[,,y] <- replace(b3[,,y], bqa[,,y] %in% bad_pixel_flag, NA)
        b4[,,y] <- replace(b4[,,y], bqa[,,y] %in% bad_pixel_flag, NA)
        b5[,,y] <- replace(b5[,,y], bqa[,,y] %in% bad_pixel_flag, NA)
        b7[,,y] <- replace(b7[,,y], bqa[,,y] %in% bad_pixel_flag, NA)
    }

    # for landsat 5-7                           # For landsat 8
    ndvi.array <- (b4 - b3)/(b4 + b3)           #ndvi.array <- (b5 - b4)/(b5 + b4)
    # for landsat 5-7
    ndbi.array <- (b5 - b4)/(b5 + b4)           # for ladnsat 8: ndbi.array <- (b6 - b5)/(b6 + b5)

    #calcualte seasonal medians
    s1_ind <- which(mth <= 3)                           # s1 includes months of Jan, Feb, March,
    s2_ind <- which(mth <= 6 & mth > 3)                 # s2 includes months of Apr, May, Jun
    s3_ind <- which(mth <= 8 & mth > 6)                 # s3 includes months of Jul, Aug, Sep
    s4_ind <- which(mth <= 12 & mth >8)                 # s4 includes months of Oct, Nov, Dec
    seasons_ind <- list(s1 = s1_ind, s2 = s2_ind, s3 = s3_ind, s4 = s4_ind)

    # get ndvi seasonal medians
    array.dim <- dim(b1)
    mat.nrow <- array.dim[1]*array.dim[2]
    mat.ncol <- array.dim[3]

    ndvi.mat <- array(ndvi.array, dim = c(mat.nrow, mat.ncol))
    ndbi.mat <- array(ndvi.array, dim = c(mat.nrow, mat.ncol))
    b1.mat <- array(b1, dim = c(mat.nrow, mat.ncol))
    b2.mat <- array(b2, dim = c(mat.nrow, mat.ncol))
    b3.mat <- array(b3, dim = c(mat.nrow, mat.ncol))
    b4.mat <- array(b4, dim = c(mat.nrow, mat.ncol))
    b5.mat <- array(b5, dim = c(mat.nrow, mat.ncol))
    b7.mat <- array(b7, dim = c(mat.nrow, mat.ncol))

    
    # get seasonal medians for every band
    allBands_matrix <- list(b1 = b1.mat, b2 = b2.mat, b3 = b3.mat, b4 = b4.mat, b5 = b5.mat, b7 = b7.mat, 
                            ndvi = ndvi.mat, ndbi = ndbi.mat)


    s1_bands <- lapply(allBands_matrix, FUN = function(t) {
        tmp <- t[, seasons_ind[['s1']]]
        return(tmp)
    })

    s1_medians <- lapply(s1_bands, FUN = function(k){
        if(length(dim(k)) == 0){
            return(k)
        }else{
            k <- rowMedians(k, na.rm = TRUE)
        }
    })
    
    s2_bands <- lapply(allBands_matrix, FUN = function(t) {
        tmp <- t[, seasons_ind[['s2']]]
        return(tmp)
    })

    s2_medians <- lapply(s2_bands, FUN = function(k){
        if(length(dim(k)) == 0){
            return(k)
        }else{
            k <- rowMedians(k, na.rm = TRUE)
        }
    })

    s3_bands <- lapply(allBands_matrix, FUN = function(t) {
        tmp <- t[, seasons_ind[['s3']]]
        return(tmp)
    })

    s3_medians <- lapply(s3_bands, FUN = function(k){
        if(length(dim(k)) == 0){
            return(k)
        }else{
            k <- rowMedians(k, na.rm = TRUE)
        }
    })
    
    s4_bands <- lapply(allBands_matrix, FUN = function(t) {
        tmp <- t[, seasons_ind[['s4']]]
        return(tmp)
    })

    s4_medians <- lapply(s4_bands, FUN = function(k){
        if(length(dim(k)) == 0){
            return(k)
        }else{
            k <- rowMedians(k, na.rm = TRUE)
        }
    })

    s4_medians_array <- array(as.numeric(unlist(s4_medians)), dim = c(nc.nrow, nc.ncol, 8))
    s1_medians_array <- array(as.numeric(unlist(s1_medians)), dim = c(nc.nrow, nc.ncol, 8))
    s2_medians_array <- array(as.numeric(unlist(s2_medians)), dim = c(nc.nrow, nc.ncol, 8))
    s3_medians_array <- array(as.numeric(unlist(s3_medians)), dim = c(nc.nrow, nc.ncol, 8))
    
    output.path  <- "D:/Zekun/Imp_classification/ard_season_medians/"
    output.name <- paste(output.path, 'LT_05_', yr, '_seasonal_medians.nc', sep = "")
    layers <- names(allBands_matrix)
    
    data_name <- 'sr'
    londim <- ncdim_def("lon","degrees_east", as.double(lon.in))
    latdim <- ncdim_def("lat", "degrees_north", as.double(lat.in))
    layerdim <- ncdim_def("layer", "layer", as.double(1:8))
    
    fillvalue <- NA
    data_long_name <- "seasonal median of surface_reflectance"
    s1_def <- ncvar_def(name = "s1_med", units = "", dim = list(londim, latdim, layerdim), fillvalue, data_long_name, prec = "single")
    s2_def <- ncvar_def(name = "s2_med", units = "", dim = list(londim, latdim, layerdim), fillvalue, data_long_name, prec = "single")
    s3_def <- ncvar_def(name = "s3_med", units = "", dim = list(londim, latdim, layerdim), fillvalue, data_long_name, prec = "single")
    s4_def <- ncvar_def(name = "s4_med", units = "", dim = list(londim, latdim, layerdim), fillvalue, data_long_name, prec = "single")
    
    ncf.out <- nc_create(output.name, list(s1_def, s2_def, s3_def, s4_def), force_v4 =TRUE)
    ncatt_put(ncf.out, "lon", "axis", "X")
    ncatt_put(ncf.out, "lat", "axis", "Y")
    ncatt_put(ncf.out, "layer", "axis", "L")
    ncatt_put(ncf.out, 0, "CRS", crs.in$value)
    ncatt_put(ncf.out, 0, "cell_size", cell.in$value)
    ncatt_put(ncf.out, 0, "xmin", xmin.in$value)
    ncatt_put(ncf.out, 0, "xmax", xmax.in$value)
    ncatt_put(ncf.out, 0, "ymin", ymin.in$value)
    ncatt_put(ncf.out, 0, "ymax", ymax.in$value)
    ncatt_put(ncf.out, 0, 'nrow ncol', "2814 by 2115")
    ncatt_put(ncf.out, 0, 'layers of each season', c("b1, b2, b3, b4, b5, b7, ndvi, ndbi"))

    ncvar_put(ncf.out, s1_def, s1_medians_array)
    ncvar_put(ncf.out, s2_def, s2_medians_array)
    ncvar_put(ncf.out, s3_def, s3_medians_array)
    ncvar_put(ncf.out, s4_def, s4_medians_array)

    nc_close(ncf.out)
}


s2.b1.ref <- b1.mat[, seasons_ind[['s2']]]
s2.b1.ref.med <- rowMedians(s2.b1.ref, na.rm = TRUE)
s2.b1.out <- s2_medians_array[,,1]
s2.b2.out <- s2_medians_array[,,2]
s2.b2.ref <- b2.mat[, seasons_ind[['s2']]]
s2.b2.ref.med <- rowMedians(s2.b2.ref, na.rm = TRUE)
diff1 <- s2.b1.out - s2.b1.ref.med
diff2 <- s2.b2.out - s2.b2.ref.med
mean(diff1, na.rm =TRUE)
mean(diff2, na.rm = TRUE)
##############################################################################################
