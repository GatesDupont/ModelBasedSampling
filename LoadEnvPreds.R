# Gates Dupont #
# F18          #
# # # # # # # #

library(raster)
library(FedData)
library(dismo)
library(dplyr)

"NOTES"
# Include elevation


#------------------------------------1. CropScape------------------------------------
crops_str = "~/Remote Sensing Data/CropScapeBLRA/CDL_2017_clip_20181016131456_984201100.tif"
crops.raw = raster(crops_str)

#----Adjusting the extent of CropScape----

crops.ext = extent(crops.raw)

# Fix ymin
crops.ext@ymin = crops.ext@ymin+((crops.ext@ymax - crops.ext@ymin)*0.1)
crops = crop(crops.raw, crops.ext)

# Fix ymax
crops.ext@xmin = crops.ext@xmin+((crops.ext@xmax - crops.ext@xmin)*0.1)
crops = crop(crops, crops.ext)


#------------------------------------2. NLCD------------------------------------
label = "CalNLCD"
if(F){
  nlcd = get_nlcd(template=crops, label=label, year = 2011, 
                  dataset = "landcover", 
                  raw.dir = "~/RAW/Remote Sensing Data/NLCD", 
                  extraction.dir = paste0("~/RAW/Remote Sensing Data/NLCD",
                                          label, "/NLCD 2011"))
}


#------------------------------------3. BLRA OCC------------------------------------

#----Download from gbif----
#occ = gbif("Laterallus", "jamaicensis*", geo=TRUE, ext=extent(crops))
pres.raw = occ[occ$adm1 == "California", c("lat", "lon")]
pres = pres.raw[complete.cases(pres.raw),]

#----Convert to sp---
coordinates(pres) = ~lon+lat
crs(pres) = CRS("+init=epsg:4326")


#----Crop occurence data----
pres = pres %>%
  spTransform(crs(crops)) %>%
  crop(crops)
plot(crops); points(pres, col="red", pch=20, cex=3)

#----Removing duplicates----
pres.nD = as.data.frame(pres)
pres.nD$combined = paste0(pres.nD$lat, ", ", pres.nD$lon)
pres.nD = pres.nD %>%
  distinct(lon, lat, .keep_all = T) # Keeping only unique values

#----Back to spatial----
pres.nD.sp = pres.nD[,c(1,2)]
coordinates(pres.nD.sp) = ~lon+lat
crs(pres.nD.sp) = crs(crops)

#----Pseudo-absence points----
set.seed(4797)
bg = as.data.frame(randomPoints(crops, 750)) # Not subsetting in polygons because small area
points(bg, col="white")

#---Combining pa----
blra = as.data.frame(pres.nD.sp)
blra$pa = 1
bg$pa = 0
colnames(bg) = c("lon", "lat", "pa")
blra = rbind(blra, bg)


#------------------------------------4. Extract CropScape------------------------------------

#------------------------------------5. Extract NLCD------------------------------------

#------------------------------------6. Model------------------------------------

#------------------------------------7. Prediction grid------------------------------------

#------------------------------------8. Model predictions------------------------------------

#------------------------------------9. Model predictions------------------------------------
