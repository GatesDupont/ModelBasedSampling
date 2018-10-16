# Gates Dupont #
# F18          #
# # # # # # # #

library(raster)
library(FedData)

#----Loading CropScape data----
crops_str = "~/Remote Sensing Data/CropScapeBLRA/CDL_2017_clip_20181016131456_984201100.tif"
crops = raster(crops_str)

#----nlcd----
label = "CalNLCD"
#nlcd = get_nlcd(template=crops, label=label, year = 2011, 
                dataset = "landcover", 
                raw.dir = "~/RAW/Remote Sensing Data/NLCD", 
                extraction.dir = paste0("~/RAW/Remote Sensing Data/NLCD",
                                        label, "/NLCD 2011"))
