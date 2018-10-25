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


#------------------------------------0. Extent------------------------------------

#----Hand-selected points from gmaps----
ur = c(41.056109, -120.771108)
ll = c(38.533879, -122.549340)

#----Spatial points object with embedded extent----
corners.coords = t(data.frame(rev(ur), rev(ll)))
study.extent.corners = SpatialPoints(corners.coords, CRS("+init=epsg:4326"))

#------------------------------------1. Get CropScape------------------------------------

#----Pulling cropscape data----
CropScape = getCDL("CA", 2017)
crops.raw = CropScape$CA2017

#----Cropping cropscape data to study extent----
crops = crop(crops.raw, spTransform(study.extent.corners, crs(CropScape$CA2017)))

#------------------------------------2. Get NLCD------------------------------------

#----Puilling nlcd data----
label = "CalNLCD"
nlcd = get_nlcd(template=crops, label=label, year = 2011) # auto-crops

#------------------------------------3. Get Elevation------------------------------------

#----Pulling elevation data----
srtm1 = getData("SRTM", lat=ur[1], lon=ur[2])
srtm2 = getData("SRTM", lat=ll[1], lon=ll[2])

#----Stitching elevation tiles----
elevation = mosaic(srtm1, srtm2, fun=mean)
elevation = crop(elevation, spTransform(study.extent.corners, crs(elevation)))

#------------------------------------4. Get WorldClim------------------------------------
climate = getData('worldclim', var='bio', res=0.5, lat=ur[1], lon=ur[2])
climate = crop(climate, spTransform(study.extent.corners, crs(climate)))

#------------------------------------5. Get eBird------------------------------------

#----Download from gbif----
occ = gbif("Laterallus", "jamaicensis*", geo=TRUE, ext=extent(study.extent.corners))
pres.raw = occ[occ$adm1 == "California", c("lat", "lon")]
pres = pres.raw[complete.cases(pres.raw),]

#----Removing duplicates----
pres.nD = as.data.frame(pres)
pres.nD$combined = paste0(pres.nD$lat, ", ", pres.nD$lon)
pres.nD = pres.nD %>%
  distinct(lon, lat, .keep_all = T) # Keeping only unique values
pres.nD = pres.nD[,c(1,2)]

#----Pseudo-absence points----
set.seed(4797)
bg = as.data.frame(randomPoints(crops, 750)) # Not subsetting in polygons because small area
colnames(bg) = c("lon", "lat")
coordinates(bg) = ~lon+lat
crs(bg) = crs(crops)
bg = spTransform(bg, CRS("+init=epsg:4326"))
bg = data.frame(bg@coords)

#----Combining presence/absence----
species = data.frame(rbind(pres.nD, bg))
species$pa = c(rep(1, length(pres.nD$lat)), rep(0, length(bg$lat)))
species.spatial = SpatialPointsDataFrame(cbind(species$lon, species$lat), data=species, proj4string = CRS("+init=epsg:4326"))

#------------------------------------6. Extract CropScape------------------------------------

#----Converting data to CropScape crs----
species.spatial.crops = spTransform(species.spatial, crs(crops))

#----Extracting CropScape values----
crops.vx = velox(stack(crops))
spol = gBuffer(species.spatial.crops, width=100, byid=TRUE)
spdf = SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
ex.mat = crops.vx$extract(spdf)
rm(crops.vx)

#----Calculating proportional cover----
pb = txtProgressBar(min = 1, max = length(ex.mat), initial = 1) 
unique.crops = sort(unique(values(crops)))
prop.rep = rep(0,74)
if(T){
  if(exists("prop.lc.df")){rm(prop.lc.df)}
  prop.lc.df = data.frame(1:74)
  for(i in 1:length(ex.mat)){
    setTxtProgressBar(pb,i)
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
colnames(prop.lc.df) = paste0("cs", sort(unique(values(crops))))
Prop.CropScape = prop.lc.df

#------------------------------------7. Extract NLCD------------------------------------

#----Converting data to nlcd crs----
species.spatial.nlcd = spTransform(species.spatial, crs(nlcd))

#----Extracting NLCD values----
nlcd.vx = velox(stack(nlcd))
spol = gBuffer(species.spatial.nlcd, width=100, byid=TRUE)
spdf = SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
ex.mat = nlcd.vx$extract(spdf)
rm(nlcd.vx)

#----Calculating proportional cover----
pb = txtProgressBar(min = 1, max = length(ex.mat), initial = 1) 
unique.nlcd = sort(unique(values(nlcd))) ###
prop.rep = rep(0,16)
if(T){
  if(exists("prop.lc.df")){rm(prop.lc.df)}
  prop.lc.df = data.frame(1:16)
  for(i in 1:length(ex.mat)){
    setTxtProgressBar(pb,i)
    if(exists("empty.pr.lc")){rm(empty.pr.lc)}
    if(exists("lc.raw")){rm(lc.raw)}
    empty.pr.lc = data.frame(Var1=unique.nlcd, prop=prop.rep)
    lc.raw = as.data.frame(table(unlist(ex.mat[[i]])))
    lc.raw$prop = lc.raw$Freq/sum(lc.raw$Freq)
    empty.pr.lc$prop[match(lc.raw$Var1, empty.pr.lc$Var1)] <- lc.raw$prop
    proportion.lc = data.frame(empty.pr.lc$prop)
    prop.lc.df = cbind(prop.lc.df, proportion.lc) # Error pops up for length.
  }
  prop.lc.df = prop.lc.df[,-1]
  prop.lc.df = t(prop.lc.df)
  rownames(prop.lc.df) = NULL
}
colnames(prop.lc.df) = paste0("nlcd", sort(unique(values(nlcd))))
Prop.NLCD = prop.lc.df

#------------------------------------8. Extract Elevation------------------------------------

#----Converting data to elevation crs----
species.spatial.elevation = spTransform(species.spatial, crs(elevation))
elevation.extVals = extract(elevation, species.spatial.elevation)

#------------------------------------9, Extract Climate------------------------------------

#----Converting data to climate crs----
species.spatial.climate = spTransform(species.spatial, crs(climate))
climate.extVals = data.frame(extract(climate, species.spatial.climate))

#----Bringing everything together----
species.df = as.data.frame(species.spatial@coords)
colnames(species.df) = c("long", "lat")
species.df = cbind(pa = species.spatial@data[,3], species.df, Prop.CropScape, Prop.NLCD, elevation.extVals, climate.extVals)
species.df = species.df[complete.cases(species.df), ]

#------------------------------------10. Model------------------------------------

#----Splitting into training and testing----
set.seed(4797) 
sample = sample.split(species.df$pa, SplitRatio = .7)
train = subset(species.df, sample == TRUE)
test  = subset(species.df, sample == FALSE)

#----Running Random Forest model----
rf = randomForest(pa ~ ., species.df, ntree=50)

#----Evaluating model performance----
evaluate(test[test$pa != 0,], test[test$pa == 0,], rf)
