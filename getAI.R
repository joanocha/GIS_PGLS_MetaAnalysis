setwd("~/aridity_index/")
library(rgdal)  
library(raster) 

ai <- raster(readGDAL("~/aridity_index/AI_annual/ai_yr/w001001.adf"))

dir <- "~/path/to/IUCN_shapefiles/" # make sure you downloaded the shapefiles for species distribution range in IUCN
sps <- list.files(dir, "*shp$", recursive = TRUE) 

data.sp <- data.frame(name=sps, min=NA, max=NA, mean=NA, sd=NA)

#loop for each species
for (i in 1:length(sps)) {    
  sp.area <- readOGR(file.path(dir, sps[1]))
  
  # extract AI
  sp.ai <- extract(ai, sp.area)
  sp.ai.vector <- as.numeric(unlist(sp.ai))
  
  data.sp[i,"min"] <- min(sp.ai.vector, na.rm=TRUE)
  data.sp[i,"max"] <- max(sp.ai.vector, na.rm=TRUE)
  data.sp[i,"mean"] <- mean(sp.ai.vector, na.rm=TRUE)
  data.sp[i,"sd"] <- sd(sp.ai.vector, na.rm=TRUE)

}

write.csv(data.sp, "AI.sp.csv") 


