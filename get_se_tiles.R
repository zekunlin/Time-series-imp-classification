

library(rgdal)
setwd("D:/Zekun/Landsat_ARD/Landsat_ARD/ARD_Uni_shapefiles/")
library(rgdal)
ard_shp <- readOGR("D:/Zekun/Dissertation/Data/se_us_shp/se_ard_tiles.shp")

for (i in 1:length(ard_shp)) {
    plgon_i <- ard_shp[i, ]
    shp_name <- paste("h", plgon_i@data$h, "_", "v", plgon_i@data$v, sep = "")
    writeOGR(
        obj = plgon_i,
        dsn = paste("D:/Zekun/Landsat_ARD/Landsat_ARD/ARD_Uni_shapefiles/", shp_name, ".shp", sep = ""),
        driver = "ESRI Shapefile", layer = shp_name
    )
    zip.list <- list.files("D:/Zekun/Landsat_ARD/Landsat_ARD/ARD_Uni_shapefiles/",
        pattern = shp_name, full.names = T
    )

    zip(zipfile = paste(shp_name, ".zip", sep = ""), files = zip.list)
}