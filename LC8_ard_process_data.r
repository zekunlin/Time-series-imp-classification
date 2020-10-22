getwd()
library(raster)
library(ncdf4)
library(foreach)
library(stringr)
library(doParallel)

data_repo <- "D:/Zekun/Landsat_ARD/Landsat_ARD/wake_ard_period3/"

## explore and store data in netcdf format

b2 <- list.files(data_repo, pattern = "^.*SRB2.*", full.names = T)
b3 <- list.files(data_repo, pattern = "^.*SRB3.*", full.names = T)
b4 <- list.files(data_repo, pattern = "^.*SRB4.*", full.names = T)
b5 <- list.files(data_repo, pattern = "^.*SRB5.*", full.names = T)
b6 <- list.files(data_repo, pattern = "^.*SRB6.*", full.names = T)
b7 <- list.files(data_repo, pattern = "^.*SRB7.*", full.names = T)
bqa <- list.files(data_repo, pattern = "^.*PIXELQA.doy*", full.names = T)


for(yr in 2013:2018) {
    
    
    b2 <- list.files(data_repo, pattern = paste("^.*SRB2.doy", yr, sep = ""), full.names = TRUE)
    b3 <- list.files(data_repo, pattern = paste("^.*SRB3.doy", yr, sep = ""), full.names = TRUE)
    b4 <- list.files(data_repo, pattern = paste("^.*SRB4.doy", yr, sep = ""), full.names = TRUE)
    b5 <- list.files(data_repo, pattern = paste("^.*SRB5.doy", yr, sep = ""), full.names = TRUE)
    b6 <- list.files(data_repo, pattern = paste("^.*SRB6.doy", yr, sep = ""), full.names = TRUE)
    b7 <- list.files(data_repo, pattern = paste("^.*SRB7.doy", yr, sep = ""), full.names = TRUE)
    bqa <- list.files(data_repo, pattern = paste("^.*PIXELQA.doy", yr, sep = ""), full.names = TRUE)

    b2_stack <- stack(b2)
    b2.array <- array(b2_stack, dim = dim(b2_stack))
    b3_stack <- stack(b3)
    b3.array <- array(b3_stack, dim = dim(b3_stack))
    b4_stack <- stack(b4)
    b4.array <- array(b4_stack, dim = dim(b4_stack))
    b5_stack <- stack(b5)
    b5.array <- array(b5_stack, dim = dim(b5_stack))
    b6_stack <- stack(b6)
    b6.array <- array(b6_stack, dim = dim(b6_stack))
    b7_stack <- stack(b7)
    b7.array <- array(b7_stack, dim = dim(b7_stack))
    qa_stack <- stack(bqa)
    qa.array <- array(qa_stack, dim = dim(qa_stack))

    
    lon <- xFromCol(b2_stack)
    lat <- yFromRow(b2_stack)
    time <- str_extract(b2, '[0-9][0-9][0-9][0-9][0-9][0-9][0-9]')
    crs.info <- as.character(crs(b2_stack))
    extent.info <-  extent(b2_stack)
    res.info <- as.character(res(b2_stack))
    varname <- 'sr'
    Date.range <- c(time[1], tail(time, n=1))

    ncpath <- "D:/Zekun/Imp_classification/ard_nc_files/"
    ncname <- paste("LC08_Wake_ncdf4_", yr, ".nc", sep = "")
    file_name <- paste(ncpath, ncname, sep = "")
    data_name <- 'sr'
    londim <- ncdim_def("lon","degrees_east", as.double(lon))
    latdim <- ncdim_def("lat", "degrees_north", as.double(lat))
    timedim <- ncdim_def("time", "day_of_year", as.double(time))
    fillvalue <- NA
    data_long_name <- "surface_reflectance"

    b2_def <- ncvar_def(name = "b2_sr", units = "", dim = list(londim, latdim, timedim), fillvalue, data_long_name, prec = "single")
    b3_def <- ncvar_def(name = "b3_sr", units = "", dim = list(londim, latdim, timedim), fillvalue, data_long_name, prec = "single")
    b4_def <- ncvar_def(name = "b4_sr", units = "", dim = list(londim, latdim, timedim), fillvalue, data_long_name, prec = "single")
    b5_def <- ncvar_def(name = "b5_sr", units = "", dim = list(londim, latdim, timedim), fillvalue, data_long_name, prec = "single")
    b6_def <- ncvar_def(name = "b6_sr", units = "", dim = list(londim, latdim, timedim), fillvalue, data_long_name, prec = "single")
    b7_def <- ncvar_def(name = "b7_sr", units = "", dim = list(londim, latdim, timedim), fillvalue, data_long_name, prec = "single")
    bqa_def <- ncvar_def(name = "qa", units = "", dim = list(londim, latdim, timedim), fillvalue, data_long_name, prec = "single")

    ncout <- nc_create(file_name, list(b2_def, b3_def, b4_def, b5_def, b6_def, b7_def, bqa_def), force_v4 =TRUE)
    ncatt_put(ncout, "lon", "axis", "X")
    ncatt_put(ncout, "lat", "axis", "Y")
    ncatt_put(ncout, "time", "axis", "T")
    ncatt_put(ncout, 0, "CRS", crs.info)
    ncatt_put(ncout, 0, "cell_size", res.info)
    ncatt_put(ncout, 0, "xmin", extent.info[1])
    ncatt_put(ncout, 0, "xmax", extent.info[2])
    ncatt_put(ncout, 0, "ymin", extent.info[3])
    ncatt_put(ncout, 0, "ymax", extent.info[4])

    
    ncvar_put(ncout, b2_def, b2.array)
    ncvar_put(ncout, b3_def, b3.array)
    ncvar_put(ncout, b4_def, b4.array)
    ncvar_put(ncout, b5_def, b5.array)
    ncvar_put(ncout, b6_def, b6.array)
    ncvar_put(ncout, b7_def, b7.array)
    ncvar_put(ncout, bqa_def, qa.array)

    nc_close(ncout)
}

ncf <- nc_open(file_name)
nc_close(ncf)

lat_import <- ncvar_get(ncf, "lat")
lon_import <- ncvar_get(ncf, "lon")
sr_import <- ncvar_get(ncf, "sr")
time_import <- ncvar_get(ncf, "time")
ncraster <- brick(file_name, varname = "b6_sr")
#time <- str_extract(b1, '[0-9][0-9][0-9][0-9][0-9][0-9][0-9]') # get doy of each image