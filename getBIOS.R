setwd("~/WorldClim/")
library(rgdal)  
library(raster)
bio<-raster("~/WorldClim/wc2.1_30s_bio/wc2.1_30s_bio_1.tif") # example for BIO1, you can set as many variables as Temperatures
dir <- "~/path/to/IUCN_shapefiles/" # make sure you downloaded the shapefiles for species distribution range in IUCN
sps <- list.files(dir, "*shp$", recursive = TRUE) 

data.sp <- data.frame(name=sps, min=NA, max=NA, mean=NA, sd=NA)

#loop per species
for (i in 1:length(sps)) {    
  # ler o vector
  sp.area <- readOGR(file.path(dir, sps[i]))
  
  # extract BIO data per species distribution reange
  sp.bio <- extract(bio, sp.area)
  sp.bio.vector <- as.numeric(unlist(sp.bio))
  
  
  data.sp[i,"min"] <- min(sp.bio.vector, na.rm=TRUE)
  data.sp[i,"max"] <- max(sp.bio.vector, na.rm=TRUE)
  data.sp[i,"mean"] <- mean(sp.bio.vector, na.rm=TRUE)
  data.sp[i,"sd"] <- sd(sp.bio.vector, na.rm=TRUE)
}

write.csv(data.sp, "bio1.csv") 
