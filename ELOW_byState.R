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
library(leaflet)
library(maps)
library(viridis)


#------------------------------------0. Extent------------------------------------

#----State Contour----
states.full = c("Texas")
us = getData('GADM', country = 'US', level = 1)
st.contour = us[us$NAME_1 %in% states.full,]
st.contour = spTransform(st.contour, CRS("+init=epsg:4326"))
st.contour = aggregate(st.contour) # remove boundaries between states

#----Hand-selected points from gmaps----
ur = c(bbox(st.contour)[2,2], bbox(st.contour)[1,2])
ll = c(bbox(st.contour)[2,1], bbox(st.contour)[1,1])
ul = c(ur[1],ll[2])
lr = c(ll[1],ur[2])

extent.long = c(ur[2],ll[2],ul[2],lr[2])
extent.lat = c(ur[1],ll[1],ul[1],lr[1])
extent.df = data.frame(cbind(extent.long, extent.lat))

#----Spatial points object with embedded extent----
coordinates(extent.df) = ~extent.long+extent.lat
study.extent.corners = SpatialPoints(extent.df, CRS("+init=epsg:4326"))

#------------------------------------1. Get CropScape------------------------------------

#----Pulling cropscape data----
CropScape = getCDL("TX", 2017)
crops.raw = CropScape$TX2017

#----Cropping cropscape data to study extent----
crops = crop(crops.raw, extent(spTransform(study.extent.corners, crs(crops.raw)))*1.1)

#------------------------------------2. Get NLCD------------------------------------

#----Puilling nlcd data----
label = "CalNLCD"
nlcd = get_nlcd(template=crops, label=label, year = 2011, force.redo = T) # auto-crops

#------------------------------------3. Get Elevation------------------------------------

#----Pulling elevation data----
lat.increment = (ur[1]-ll[1])/2
long.increment = (ur[2]-ll[2])/3
e.lats = c(ul[1], ul[1]-lat.increment, ul[1]-(2*lat.increment))
e.longs = c(ul[2], ul[2]+long.increment, ul[2]+(2*long.increment), ul[2]+(3*long.increment))

srtm1 = getData("SRTM", lat=e.lats[1], lon=e.longs[1])
srtm2 = getData("SRTM", lat=e.lats[1], lon=e.longs[2])
srtm3 = getData("SRTM", lat=e.lats[1], lon=e.longs[3])
srtm4 = getData("SRTM", lat=e.lats[1], lon=e.longs[4])
srtm5 = getData("SRTM", lat=e.lats[2], lon=e.longs[1])
srtm6 = getData("SRTM", lat=e.lats[2], lon=e.longs[2])
srtm7 = getData("SRTM", lat=e.lats[2], lon=e.longs[3])
srtm8 = getData("SRTM", lat=e.lats[2], lon=e.longs[4])
srtm9 = getData("SRTM", lat=e.lats[3], lon=e.longs[1])
srtm10 = getData("SRTM", lat=e.lats[3], lon=e.longs[2])
srtm11 = getData("SRTM", lat=e.lats[3], lon=e.longs[3])
srtm12 = getData("SRTM", lat=e.lats[3], lon=e.longs[4])

#----Stitching elevation tiles----
elevation = mosaic(srtm1, srtm2, srtm3, srtm4, srtm5, srtm6, srtm7, srtm8, srtm9, srtm10, srtm11, srtm12, fun=max)
rm(srtm1, srtm2, srtm3, srtm4, srtm5, srtm6, srtm7, srtm8, srtm9, srtm10, srtm11, srtm12)
elevation = crop(elevation, extent(spTransform(study.extent.corners, crs(elevation)))*1.1)

#------------------------------------4. Get WorldClim------------------------------------
climate1 = getData('worldclim', var='bio', res=0.5, lat=ll[1], lon=ll[2])
climate2 = getData('worldclim', var='bio', res=0.5, lat=ur[1], lon=ur[2])
climate = mosaic(climate1, climate2, fun=mean)

climate = crop(climate, extent(spTransform(study.extent.corners, crs(climate)))*1.1)

#------------------------------------5. Get eBird------------------------------------

#----Download from gbif----
occ = gbif("Micrathene", "whitneyi*", geo=TRUE, ext=extent(study.extent.corners))
pres.raw = occ[occ$adm1 == "Texas", c("lat", "lon")]
pres = pres.raw[complete.cases(pres.raw),]

#----Removing duplicates----
pres.nD = as.data.frame(pres)
pres.nD$combined = paste0(pres.nD$lat, ", ", pres.nD$lon)
pres.nD = pres.nD %>%
  distinct(lon, lat, .keep_all = T) # Keeping only unique values
pres.nD = pres.nD[,c(1,2)]

#----Grid Sampling Presence---
r = raster(study.extent.corners)
res(r) = 0.075
r = extend(r, extent(r)+1)

coordinates(pres.nD) <- ~lon+lat
pres.samp = gridSample(pres.nD, r, n=1)
points(pres.samp, pch=20, col="purple", cex=0.3)

p = rasterToPolygons(r)
plot(p)
points(pres.nD$lat~pres.nD$lon, pch=20, col="red", cex=0.3)

pres.nD = as.data.frame(pres.samp)

#----Pseudo-absence points----
set.seed(4797)
bg = as.data.frame(spsample(st.contour,n=length(pres.nD$lat)*10,"random"))
colnames(bg) = c("lon", "lat")
coordinates(bg) = ~lon+lat
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
spol = gBuffer(species.spatial.crops, width=1000, byid=TRUE)
spdf = SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
ex.mat = crops.vx$extract(spdf)
rm(crops.vx)

#----Calculating proportional cover----
pb = txtProgressBar(min = 1, max = length(ex.mat), initial = 1) 
unique.crops = sort(unique(values(crops)))
prop.rep = rep(0,length(unique.crops))
if(T){
  if(exists("prop.lc.df")){rm(prop.lc.df)}
  prop.lc.df = data.frame(1:length(unique.crops))
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
colnames(prop.lc.df) = paste0("cs", unique.crops)
Prop.CropScape = prop.lc.df

#------------------------------------7. Extract NLCD------------------------------------

#----Converting data to nlcd crs----
species.spatial.nlcd = spTransform(species.spatial, crs(nlcd))

#----Extracting NLCD values----
nlcd.vx = velox(stack(nlcd))
spol = gBuffer(species.spatial.nlcd, width=1000, byid=TRUE)
spdf = SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
ex.mat = nlcd.vx$extract(spdf)
rm(nlcd.vx)

#----Calculating proportional cover----
pb = txtProgressBar(min = 1, max = length(ex.mat), initial = 1) 
unique.nlcd = sort(unique(values(nlcd)))
prop.rep = rep(0,length(unique.nlcd))
if(T){
  if(exists("prop.lc.df")){rm(prop.lc.df)}
  prop.lc.df = data.frame(1:length(unique.nlcd))
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
colnames(prop.lc.df) = paste0("nlcd", unique.nlcd)
Prop.NLCD = prop.lc.df

#------------------------------------8. Extract Elevation------------------------------------

#----Converting data to elevation crs----
species.spatial.elevation = spTransform(species.spatial, crs(elevation))
elevation.extVals = extract(elevation, species.spatial.elevation)

#------------------------------------9. Extract Climate------------------------------------

#----Converting data to climate crs----
species.spatial.climate = spTransform(species.spatial, crs(climate))
climate.extVals = data.frame(extract(climate, species.spatial.climate))

#----Bringing everything together----
species.df = as.data.frame(species.spatial@coords)
colnames(species.df) = c("long", "lat")
species.df = cbind(pa = species.spatial@data[,3], species.df, Prop.CropScape, Prop.NLCD, elevation.extVals, climate.extVals)
species.df = species.df[complete.cases(species.df), ]

#------------------------------------10. Model------------------------------------

#----Making 10 data sets----
rf.df.pres = species.df[species.df$pa == 1,]
rf.dfs = list()
starts = seq(length(rf.df.pres$pa)+1,length(species.df$pa), by=length(rf.df.pres$pa))
for(i in 1:10){
  rf.dfs[[i]] =  data.frame(rbind(rf.df.pres, species.df[c(starts[i]:(starts[i]+length(rf.df.pres$pa))),]))
  rf.dfs[[i]]= rf.dfs[[i]][1:(length(rf.dfs[[i]]$pa)-1),]
} # last one in rf.dfs[[10]] is NA, just a simple counting problem.

#----Loop to create training & testing sets----

#----Running Random Forest models---- # Fix to run on training
rf = vector("list", 10)
for(i in 1:10){
  rf[[i]] = randomForest(pa ~ ., rf.dfs[[i]])
}

#----Loop to evaluate each rf, storing AUC----

#----AUC boxplot----

#----Splitting into training and testing----
#set.seed(4797) 
#sample = sample.split(species.df$pa, SplitRatio = .7)
#train = subset(species.df, sample == TRUE)
#test  = subset(species.df, sample == FALSE)

#----Evaluating model performance----
#evaluate(test[test$pa != 0,], test[test$pa == 0,], rf)

#------------------------------------11. Grid Coordinates------------------------------------

#----Generating prediction grid----
grid = makegrid(st.contour, cellsize = 0.05) # 0.05
grid = SpatialPoints(grid, proj4string = CRS(proj4string(st.contour)))
date() ; grid <- grid[st.contour, ] ; date()


#------------------------------------12. Grid CropScape------------------------------------

#----Converting data to CropScape crs----
species.spatial.crops = spTransform(grid, crs(crops))

#----Extracting CropScape values----
crops.vx = velox(stack(crops))
spol = gBuffer(species.spatial.crops, width=1000, byid=TRUE)
spdf = SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
ex.mat = crops.vx$extract(spdf)
rm(crops.vx, spol, spdf)

#----Calculating proportional cover----
pb = txtProgressBar(min = 1, max = length(ex.mat), initial = 1) 
unique.crops = sort(unique(crops))
prop.rep = rep(0,length(unique.crops))
if(T){
  if(exists("prop.lc.df")){rm(prop.lc.df)}
  prop.lc.df = data.frame(1:length(unique.crops))
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

colnames(prop.lc.df) = paste0("cs", unique.crops)
Prop.CropScape = prop.lc.df
rm(ex.mat, pb, unique.crops, prop.lc.df)

#------------------------------------13. Grid NLCD------------------------------------

#----Converting data to nlcd crs----
species.spatial.nlcd = spTransform(grid, crs(nlcd))

#----Extracting NLCD values----
nlcd.vx = velox(stack(nlcd))
spol = gBuffer(species.spatial.nlcd, width=1000, byid=TRUE)
spdf = SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
ex.mat = nlcd.vx$extract(spdf)
rm(nlcd.vx, spol, spdf)

#----Calculating proportional cover----
pb = txtProgressBar(min = 1, max = length(ex.mat), initial = 1) 
unique.nlcd = sort(unique(nlcd))
prop.rep = rep(0,length(unique.nlcd))
if(T){
  if(exists("prop.lc.df")){rm(prop.lc.df)}
  prop.lc.df = data.frame(1:length(unique.nlcd))
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
colnames(prop.lc.df) = paste0("nlcd", unique.nlcd)
Prop.NLCD = prop.lc.df

#------------------------------------14. Grid Elevation------------------------------------

species.spatial.elevation = spTransform(grid, crs(elevation))
elevation.extVals = extract(elevation, species.spatial.elevation)

#------------------------------------15. Grid Climate------------------------------------

#----Converting data to climate crs----
species.spatial.climate = spTransform(grid, crs(climate))
climate.extVals = data.frame(extract(climate, species.spatial.climate))

#------------------------------------16. Combining data------------------------------------
pred.df = as.data.frame(grid@coords)
colnames(pred.df) = c("long", "lat")
pred.df = cbind(pred.df, Prop.CropScape, Prop.NLCD, elevation.extVals, climate.extVals)
pred.df = pred.df[complete.cases(pred.df), ]

#------------------------------------17. Model Predictions------------------------------------

#----Model-averaged predictions----
rf.p = data.frame(predict(rf[[1]], pred.df))
for(i in 2:9){
  rf.p=cbind(rf.p, predict(rf[[i]], pred.df))
}
rf.avg = rowMeans(rf.p)

#rf.predictions = predict(rf[[1]], pred.df)
predictions = pred.df[,c("long", "lat")]
predictions$rf = rf.avg

#------------------------------------18. Plotting------------------------------------

#----Creating the raster----
SDM.raster = rasterFromXYZ(predictions)
crs(SDM.raster) = crs(grid)

#----Basic plot----
plot(rasterFromXYZ(predictions), main = "Species Distribution Model",
     col=colorRampPalette(c("darkgreen", "yellow", "red"))(100), zlim=c(0,1))

plot(rasterFromXYZ(predictions), main = "Species Distribution Model",
     col=plasma(100), zlim=c(0,1))

#----Leaflet----
pal1 = colorNumeric(rev(c("#FF0000", "#FFFF00", "#228B22")), values(SDM.raster),
                    na.color = "transparent")

pal2 = colorNumeric(colorRampPalette(c("gray", "yellow", "orange", "darkorchid", "darkorchid4", "black"))(100), values(SDM.raster),
                    na.color = "transparent")

leaflet() %>% addTiles(urlTemplate = "https://mts1.google.com/vt/lyrs=s&hl=en&src=app&x={x}&y={y}&z={z}&s=G", attribution = 'Google') %>%
  addRasterImage(SDM.raster, colors = pal1, opacity = 0.4) %>%
  addLegend(pal = pal1, values = seq(0,0.96,0.1), title = "Pr(Occurence)") %>%
  addCircleMarkers(lng=pres.nD$lon, lat = pres.nD$lat, radius=0.4)


leaflet() %>% addTiles() %>%
  addRasterImage(SDM.raster, colors = pal1, opacity = 0.5) %>%
  addLegend(pal = pal1, values = seq(0,0.96,0.1), title = "Pr(Occurence)")

ggplot() +  
  theme_void() +
  scale_fill_viridis_c() +
  xlim(min(predictions$long)-3, max(predictions$long)+3) +
  ylim(ylim=c(min(predictions$lat)-3, max(predictions$lat)+3)) +
  borders("world") +
  geom_tile(data=predictions, aes(x=predictions$long, y=predictions$lat, fill=predictions$rf), alpha=0.8)

# Implement foreach for 6 cores.
