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
colnames(prop.lc.df) = paste(sort(unique(values(crops))))
Prop.CropScape = prop.lc.df
