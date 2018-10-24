# Gates Dupont #
# F18          #
# # # # # # # #

library(raster)
library(FedData)
library(cdlTools)
library(dismo)
library(dplyr)
library(velox)
library(rgeos)
library(latticeExtra)
library(randomForest)
library(caTools)

"NOTES"
# Predictors
# https://www.gis-blog.com/r-raster-data-acquisition/

#------------------------------------0. Extent------------------------------------

#----Hand-selected points from gmaps----
ur = c(41.056109, -120.771108)
ll = c(38.533879, -122.549340)

#----Spatial points object with embedded extent----
corners.coords = t(data.frame(rev(ur), rev(ll)))
study.extent.corners = SpatialPoints(corners.coords, CRS("+init=epsg:4326"))

#------------------------------------1. CropScape------------------------------------

#----Pulling cropscape data----
CropScape = getCDL("CA", 2017)
crops.raw = raster(CropScape$CA2017)

#----Cropping cropscape data to study extent----
crops = crop(crops.raw, spTransform(study.extent.corners, crs(CropScape$CA2017)))

#------------------------------------2. NLCD------------------------------------

#----Puilling nlcd data----
label = "CalNLCD"
nlcd = get_nlcd(template=crops, label=label, year = 2011) # auto-crops

#------------------------------------3. Elevation------------------------------------

#----Pulling elevation data----
srtm1 = getData("SRTM", lat=ur[1], lon=ur[2])
srtm2 = getData("SRTM", lat=ll[1], lon=ll[2])

#----Stitching elevation tiles----
elevation = mosaic(srtm1, srtm2, fun=mean)
elevation = crop(elevation, spTransform(study.extent.corners, crs(elevation)))

#------------------------------------4. WorldClim------------------------------------
climate = getData('worldclim', var='bio', res=0.5, lat=ur[1], lon=ur[2])
climate = crop(climate, spTransform(study.extent.corners, crs(climate)))

#------------------------------------5. BLRA------------------------------------

#----Download from gbif----
occ = gbif("Laterallus", "jamaicensis*", geo=TRUE, ext=extent(study.extent.corners))
pres.raw = occ[occ$adm1 == "California", c("lat", "lon")]
pres = pres.raw[complete.cases(pres.raw),]

#----Removing duplicates----
pres.nD = as.data.frame(pres)
pres.nD$combined = paste0(pres.nD$lat, ", ", pres.nD$lon)
pres.nD = pres.nD %>%
  distinct(lon, lat, .keep_all = T) # Keeping only unique values
coordinates(pres.nD) = ~lon+lat
crs(pres.nD) = CRS("+init=epsg:4326")

#----Pseudo-absence points----
set.seed(4797)
bg = as.data.frame(randomPoints(crops, 750)) # Not subsetting in polygons because small area
coordinates(bg) = ~x+y
crs(bg) = crs(crops)
bg = spTransform(bg, crs(pres.nD))
