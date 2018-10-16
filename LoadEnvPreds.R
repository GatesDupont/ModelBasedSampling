# Gates Dupont #
# F18          #
# # # # # # # #

library(raster)
library(FedData)
library(dismo)


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
blra = occ[occ$adm1 == "California", c("lat", "lon")]
blra = blra[complete.cases(blra),]

#----Pseudo-absence points----

#----Convert to sp---
coordinates(blra) = ~lon+lat
crs(blra) = CRS("+init=epsg:4326")

#----Crop occurence data----
blra = spTransform(blra, crs(crops))
blra.cropped = crop(blra, crops)
plot(crops); points(blra.cropped, col="red", pch=20, cex=3)


#------------------------------------4. Extract CropScape------------------------------------

#------------------------------------5. Extract NLCD------------------------------------

#------------------------------------6. Model------------------------------------

#------------------------------------7. Prediction grid------------------------------------

#------------------------------------8. Model predictions------------------------------------

#------------------------------------9. Model predictions------------------------------------
