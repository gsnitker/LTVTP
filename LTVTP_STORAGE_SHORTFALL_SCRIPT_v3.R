## LTVTP Storage Shortfall Script  v.3
## Grant Snitker 
## School of Human Evolution and Social Change 
## Arizona State University
## Date: 05/12/2016

### Note: Replace data inputs (in CAPITAL LETTERS) with data for specified analytical region

######## BEGIN PREAMBLE ########
# Set the working directory.
setwd("PATH TO WORKING DIRECTORY")
# Load required packages
library("parallel")
library("raster")
library("plyr")
library("utils")
last <- function(x) { return( x[length(x)] ) } # MISC function
numcores = detectCores()
######## END PREAMBLE ########

######## BEGIN PARAMETERS ########
# Force Raster to load large rasters into memory
rasterOptions(chunksize=2e+07,maxmemory=2e+08)
# Import precipitation map for analytical region (the LTVTP project uses precip. data from PRISM Climate Group, Oregon State University, http://prism.oregonstate.edu )
prcp <- brick("RAW PRCP DATA") # rasterbrick of reconstructed precipitation values from the study region and temporal period of interest
prcp.mean = brick("RAW PRCP DATA MEAN") # rasterbrick of the mean reconstructed precipitation values from the study region and temporal period of interest
prcp.sd = brick("RAW PRCP DATA SD") # rasterbrick of the standard deviation of reconstructed precipitation values from the study region and temporal period of interest

year.range =  (900:1000)# this is the year range of interest from the precip reconstruction
year.range.buffer = (year.range[1] - (3)):(last(year.range)) # provides a buffer to calculate moving averages for the years preceding the first year in the year.range
moving.sum.count <-count(year.range[1])
######## END PARAMETERS ########



######## FUNCTION 1 #######################
###########################################
######## MOVING AVG FUCTION ########
# Append 3 years to the beginning of the laoded prcp.mean and prcp.sd
options(warn=-1) # suppress warnings for the following expression
add.range = (((year.range)-3):((year.range)-1))
options(warn=0) # reinstate warnings 
moving.value = 60 # total number of years for moving average window (in this case, 60 indicates 30 years before and 30 years after focal year)

#moving avg function
moving.avg.function = function(add.range){
  moving.avg.range = ((add.range - moving.value):add.range)
  prcp.sub = subset(prcp,moving.avg.range)
  prcp.sub.mean <- calc(prcp.sub,mean)
  prcp.sub.sd <- calc(prcp.sub, sd)
  return(list(prcp.sub.mean, prcp.sub.sd))
}
#run moving avg function for both the mean and sd maps
ptm <- proc.time()
moving.avg.nested = mclapply(add.range, function(i){moving.avg.function(i)}, mc.cores = numcores)
proc.time() - ptm

#process reults by splitting aprt mean and sd maps and aving tehm as individualbricks for the year.range
moving.avg.raster = unlist(moving.avg.nested)
prcp.mean.add<- brick(moving.avg.raster[1:length(moving.avg.raster) %% 2 == 1])
prcp.sd.add <- brick(moving.avg.raster[1:length(moving.avg.raster) %% 2 == 0])
names(prcp.mean.add) = add.range
names(prcp.sd.add) = add.range

# create a dummy raster to add to result for future functions
options(warn=-1)
dummy.raster = subset(prcp,(1:(year.range-4)))
options(warn=0)
dummy.raster[] = NA
prcp.mean.storage = addLayer(dummy.raster, prcp.mean.add, prcp.mean)
prcp.sd.storage = addLayer(dummy.raster, prcp.sd.add, prcp.sd)

######## FUNCTION 1 #######################
###########################################
######## THRESHOLD BINARY FUNCTION ########
theshold.binary.function = function(year.range.buffer){
# Moving avg and sd
  prcp.mean.storage.sub = subset(prcp.mean.storage,year.range.buffer)
  prcp.sd.storage.sub = subset(prcp.sd.storage,year.range.buffer)
# Threshold
  prcp.threshold = (prcp.mean.storage.sub - (prcp.sd.storage.sub * 0.5))
# Create below thresold maps
  prcp.instant = subset(prcp,year.range.buffer)
  below.prcp.threshold.raw = prcp.instant - prcp.threshold
  below.prcp.threshold <- reclassify(below.prcp.threshold.raw, c(-Inf,-0.00000000001,1, 0,Inf,0))  
}
######## END FUNCTION ########

######## BEGIN RUN ########
ptm <- proc.time()
below.prcp.threshold = mclapply(year.range.buffer, function(i){theshold.binary.function(i)}, mc.cores = numcores)
proc.time() - ptm
######## END RUN ########

######## BEGIN POST-PROCESSING ########
below.prcp.threshold.brick = brick(unlist(below.prcp.threshold))
names(below.prcp.threshold.brick) = year.range.buffer
below.prcp.threshold.brick = addLayer(dummy.raster, below.prcp.threshold.brick)
######## END POST-PROCESSING ########
#####################################
#####################################


######## FUNCTION 2 #######################
###########################################
######## STORAGE SHORTFALL FUNCTION ########
### New Parameters
storage.shortfall.function = function(year.range){
  # Moving sum
  moving.sum.range = (year.range - 3):(year.range) # 3 is an adjsutable value based on moving.avg.value
  below.prcp.threshold.brick.sub = subset(below.prcp.threshold.brick, moving.sum.range)
  storage.shortfall.raw <- sum(below.prcp.threshold.brick.sub)
  storage.shortfall.reclass <- reclassify(storage.shortfall.raw, c(-Inf,2,0, 2.1,Inf,1))  
  return(list(storage.shortfall.reclass))
}
######## END FUNCTION ########

######## BEGIN RUN ########
ptm <- proc.time()
storage.shortfall = mclapply(year.range, function(i){storage.shortfall.function(i)}, mc.cores = 4)
proc.time() - ptm
######## END RUN ########

######## BEGIN POST-PROCESSING ########
storage.shortfall.brick = brick(unlist(storage.shortfall))
storage.shortfall.brick
storage.shortfall.sum = sum(storage.shortfall.brick) 
count.total = count(year.range)
storage.shortfall.per = (storage.shortfall.sum) / (sum(count.total[,2]))
######## END POST-PROCESSING ########
plot(storage.shortfall.per)
######## BEGIN EXPORT ########
writeRaster(storage.shortfall.per,filename="OUTPUT PATH", format="GTiff", overwrite=TRUE)
######## END EXPORT ########
