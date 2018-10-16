# Gates Dupont #
# F18          #
# # # # # # # #

library(raster)
library(FedData)
library(dismo)
library(dplyr)
library(velox)
library(rgeos)

"NOTES"

#https://www.gis-blog.com/r-raster-data-acquisition/
# Other predictors
  # elevation
  # world clim


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


#------------------------------------3. BLRA------------------------------------

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

#----Back to spatial----
coords = cbind(blra$lon, blra$lat)
blra.coords = SpatialPoints(coords)
blra = SpatialPointsDataFrame(coords=blra.coords, data=blra)
crs(blra) = crs(crops)

#------------------------------------4. Extract CropScape------------------------------------

#----Extracting CropScape values----
crops.vx = velox(stack(crops))
spol = gBuffer(blra, width=100, byid=TRUE)
spdf = SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
ex.mat = crops.vx$extract(spdf)

#----Calculating proportional cover----
date()
if(T){
  if(exists("prop.lc.df")){rm(prop.lc.df)}
  prop.lc.df = data.frame(1:74)
  for(i in 1:length(ex.mat)){
    if(exists("empty.pr.lc")){rm(empty.pr.lc)}
    if(exists("lc.raw")){rm(lc.raw)}
    empty.pr.lc = data.frame(Var1=sort(unique(values(crops))), prop=rep(0,74))
    lc.raw = as.data.frame(table(unlist(ex.mat[[i]])))
    lc.raw$prop = lc.raw$Freq/sum(lc.raw$Freq)
    empty.pr.lc$prop[match(lc.raw$Var1, empty.pr.lc$Var1)] <- lc.raw$prop
    proportion.lc = data.frame(empty.pr.lc$prop)
    prop.lc.df = cbind(prop.lc.df, proportion.lc)
  }
  prop.lc.df = prop.lc.df[,-1]
  prop.lc.df = t(prop.lc.df)
  rownames(prop.lc.df) = NULL
}
date()
View(prop.lc.df)


#----Retransform to nlcd crs----

#------------------------------------5. Extract NLCD------------------------------------

#------------------------------------6. Model------------------------------------

#------------------------------------7. Prediction grid------------------------------------

#------------------------------------8. Model predictions------------------------------------

#------------------------------------9. Plot predictions------------------------------------
