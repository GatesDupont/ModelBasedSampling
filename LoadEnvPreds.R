# Gates Dupont #
# F18          #
# # # # # # # #

library(raster)
library(FedData)
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


#------------------------------------2. BLRA------------------------------------

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

#------------------------------------3. NLCD------------------------------------
label = "CalNLCD"
if(F){
  nlcd = get_nlcd(template=crops, label=label, year = 2011, 
                  dataset = "landcover", 
                  raw.dir = "~/RAW/Remote Sensing Data/NLCD", 
                  extraction.dir = paste0("~/RAW/Remote Sensing Data/NLCD",
                                          label, "/NLCD 2011"))
}

#------------------------------------4. Elevation------------------------------------

meanLatLon = data.frame(lon=mean(bg$lon), lat=mean(bg$lat))
coordinates(meanLatLon) = ~lon+lat
crs(meanLatLon) = crs(crops)
meanLatLon = data.frame(spTransform(meanLatLon, CRS("+init=epsg:4326")))

srtm1 = getData("SRTM", lat=meanLatLon$lat, lon=meanLatLon$lon)
srtm2 = getData("SRTM", lat=40.794833, lon=-122.151466)
elevation = mosaic(srtm1, srtm2, fun=mean)

#spplot(elevation) +
#  layer(panel.points(blra.plot@coords[,1], blra.plot@coords[,2],
#                     col="green", cex=0.001), data=blra.plot)

#------------------------------------5. WorldClim------------------------------------

climate = getData('worldclim', var='bio', res=0.5, lat=meanLatLon$lat, lon=meanLatLon$lon)

# #----Making sure this worked----
# coords = cbind(blra$lon, blra$lat)
# blra.coords = SpatialPoints(coords)
# blra = SpatialPointsDataFrame(coords=blra.coords, data=blra)
# crs(blra) = crs(nlcd)
# blra.plot = spTransform(blra, crs(climate))
# spplot(climate$bio1_11) +
#   layer(panel.points(blra.plot@coords[,1], blra.plot@coords[,2],
#                      col="green", cex=0.001), data=blra.plot)


#------------------------------------6. Extract CropScape------------------------------------

#----Extracting CropScape values----
crops.vx = velox(stack(crops))
spol = gBuffer(blra, width=100, byid=TRUE)
spdf = SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
ex.mat = crops.vx$extract(spdf)

#----Calculating proportional cover----
pb = txtProgressBar(min = 1, max = length(ex.mat), initial = 1) 
unique.crops = sort(unique(values(crops))) ###
prop.rep = rep(0,74)
date()
if(T){
  if(exists("prop.lc.df")){rm(prop.lc.df)}
  prop.lc.df = data.frame(1:74)
  for(i in 1:length(ex.mat)){
    setTxtProgressBar(pb,i)
    #print(i)
    if(exists("empty.pr.lc")){rm(empty.pr.lc)}
    if(exists("lc.raw")){rm(lc.raw)}
    empty.pr.lc = data.frame(Var1=unique.crops, prop=prop.rep)
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
colnames(prop.lc.df) = paste(sort(unique(values(crops))))
Prop.CropScape = prop.lc.df

#------------------------------------7. Extract NLCD------------------------------------

#----Retransform blra to match nlcd----
blra = spTransform(blra, crs(nlcd))


#----Extracting NLCD values----
nlcd.vx = velox(stack(nlcd))
spol = gBuffer(blra, width=100, byid=TRUE)
spdf = SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
ex.mat = nlcd.vx$extract(spdf)
rm(nlcd.vx)

#----Calculating proportional cover----
pb = txtProgressBar(min = 1, max = length(ex.mat), initial = 1) 
unique.nlcd = sort(unique(values(nlcd))) ###
prop.rep = rep(0,74)
date()
if(T){
  if(exists("prop.lc.df")){rm(prop.lc.df)}
  prop.lc.df = data.frame(1:15)
  for(i in 1:length(ex.mat)){
    setTxtProgressBar(pb,i)
    #print(i)
    if(exists("empty.pr.lc")){rm(empty.pr.lc)}
    if(exists("lc.raw")){rm(lc.raw)}
    empty.pr.lc = data.frame(Var1=unique.nlcd, prop=rep(0,15))
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
colnames(prop.lc.df) = paste(sort(unique(values(nlcd))))
Prop.NLCD = prop.lc.df
colnames(Prop.NLCD) = c("cs11", "cs21", "cs22", "cs23", "cs24", "cs31", "cs41", "cs42", 
                        "cs43", "cs52", "cs71", "cs81", "cs82", "cs90", "cs95")

blra = as.data.frame(blra)[,c(1,2, 3)]
blra = cbind(blra, Prop.CropScape, Prop.NLCD)
#write.csv(blra, "~/Model-based Sampling R/blra_CropScape_NLCD.csv")

coords = cbind(blra$lon, blra$lat)
blra.coords = SpatialPoints(coords)
blra = SpatialPointsDataFrame(coords=blra.coords, data=blra)
crs(blra) = crs(nlcd)


#------------------------------------8. Extract Elevation------------------------------------

blra = spTransform(blra, crs(elevation))
elevation.extVals = extract(elevation, blra)

#------------------------------------9, Extract Climate------------------------------------

climate.extVals = data.frame(extract(climate, blra))

#----Bringing everything together
blra = as.data.frame(blra)
blra = cbind(blra, climate.extVals)
blra$elevation = elevation.extVals
blra = blra[ , !(names(blra) %in% c("coords.x1", "coords.x2"))]
#write.csv(blra, "~/Model-based Sampling R/blra_CropScape_NLCD_WorldClim_Elev.csv")
# Coords are in crs(crops)

#------------------------------------10. Model------------------------------------

#----Splitting into training and testing----
set.seed(4797) 
sample = sample.split(blra$pa, SplitRatio = .7)
train = subset(blra, sample == TRUE)
test  = subset(blra, sample == FALSE)

rf = randomForest(pa ~ ., blra, ntree=50)

evaluate(test[test$pa != 0,], test[test$pa == 0,], rf)

#------------------------------------11. Prediction grid------------------------------------

#----Generating prediction grid -- COORDINATES----
#states.full = c("California")
# load some spatial data. Administrative Boundary
#us = raster::getData('GADM', country = 'US', level = 1)
#st.contour <- us[us$NAME_1 %in% states.full,]
#st.contour = spTransform(st.contour, crs(crops))
#grid <- makegrid(st.contour, cellsize = 100)
grid <- makegrid(bbox(crops), cellsize = 100)
grid <- SpatialPoints(grid, proj4string = CRS(proj4string(crops)))
#date() ; grid <- grid[st.contour, ] ; date()
#grid = crop(grid, crops)

#----Extracting CropScape values----
crops.vx = velox(stack(crops))
spol = gBuffer(blra, width=100, byid=TRUE)
spdf = SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
ex.mat = crops.vx$extract(spdf)

#----Calculating proportional cover----
pb = txtProgressBar(min = 1, max = length(ex.mat), initial = 1) 
unique.crops = sort(unique(values(crops))) ###
prop.rep = rep(0,74)
date()
if(T){
  if(exists("prop.lc.df")){rm(prop.lc.df)}
  prop.lc.df = data.frame(1:74)
  for(i in 1:length(ex.mat)){
    setTxtProgressBar(pb,i)
    #print(i)
    if(exists("empty.pr.lc")){rm(empty.pr.lc)}
    if(exists("lc.raw")){rm(lc.raw)}
    empty.pr.lc = data.frame(Var1=unique.crops, prop=prop.rep)
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
colnames(prop.lc.df) = paste(sort(unique(values(crops))))
Prop.CropScape = prop.lc.df

#----Reprojecting from crops to nlcd----
grid = spTransform(grid, crs(nlcd))

#----Extracting NLCD values----
nlcd.vx = velox(stack(nlcd))
spol = gBuffer(grid, width=100, byid=TRUE)
spdf = SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
ex.mat = nlcd.vx$extract(spdf)
rm(nlcd.vx)

#----Calculating proportional cover----
pb = txtProgressBar(min = 1, max = length(ex.mat), initial = 1) 
unique.nlcd = sort(unique(values(nlcd)))
prop.rep = rep(0,74)
date()
if(T){
  if(exists("prop.lc.df")){rm(prop.lc.df)}
  prop.lc.df = data.frame(1:15)
  for(i in 1:length(ex.mat)){
    setTxtProgressBar(pb,i)
    #print(i)
    if(exists("empty.pr.lc")){rm(empty.pr.lc)}
    if(exists("lc.raw")){rm(lc.raw)}
    empty.pr.lc = data.frame(Var1=unique.nlcd, prop=rep(0,15))
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
colnames(prop.lc.df) = paste(sort(unique(values(nlcd))))
Prop.NLCD = prop.lc.df
colnames(Prop.NLCD) = c("cs11", "cs21", "cs22", "cs23", "cs24", "cs31", "cs41", "cs42", 
                        "cs43", "cs52", "cs71", "cs81", "cs82", "cs90", "cs95")

blra = as.data.frame(blra)[,c(1,2, 3)]
blra = cbind(blra, Prop.CropScape, Prop.NLCD)





#------------------------------------12. Model predictions------------------------------------

#------------------------------------13. Plot predictions------------------------------------
