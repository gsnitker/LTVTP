## LTVTP Wild Plant Species Richness Script  v.5
## Grant Snitker 
## School of Human Evolution and Social Change 
## Arizona State University
## Date: 05/12/2016

### Note 1: Replace data inputs (in CAPITAL LETTERS) with data for specified analytical region
### Note 2: To use this script, the user must download MaxEnt, a java based program, from https://www.cs.princeton.edu/~schapire/maxent/.  Once downloaded, place the executable java file in the working directory for this script. 


######## BEGIN PREAMBLE ########
# Set the working directory.
setwd("PATH TO WORKING DIRECTORY")
# Load required packages
library("biomod2")
library("parallel")
library("rgdal")
library("raster")
######## END PREAMBLE ########

######## BEGIN PARAMETERS ########
### reading in species models and climate data
spNames=read.csv("DIRECTORY TO SPECIES MODELS",head=T)#species list
sp.number = (1:95) # index of species models to be used in function
year.range = (900:949) # this is the year range of interest from reconstruction

# Import precipitation map for analytical region (the LTVTP project uses precip. data from PRISM Climate Group, Oregon State University, http://prism.oregonstate.edu )
MAP.full = brick("MEAN ANNUAL PRECP DATA")
MAT.full = brick("MEAN ANNUAL TEMP DATA")

MAP = mean(subset(MAP.full, year.range)) #subset MAP and MAT to the year range of interest from reconstruction
MAT = mean(subset(MAT.full, year.range)) 

substrRight=function(x,n){#a function to help find the files associated with each species
  substr(x,nchar(x)-n+1,nchar(x))}
######## END PARAMETERS ########

######## WILD PLANT SPECIES RICHNESS FUNCTION  #######################
wild.plants.function = function(sp.number){
setwd("PATH TO WORKING DIRECTORY")#wherever your climate data and MAxEnt java file are stored
sp.i= sp.number#selecting the first species in the list
modFiles=list.files(paste(paste("./wild_plant_models/",paste(spNames[sp.i,1],spNames[sp.i,2],sep="_"),sep=""),"biomod_output/model.output",sep="/"))#the first string here is the folder that contains all of the individual species' models
modOut=modFiles[substrRight(modFiles,4)==".out"]#grabbing the model files
setwd(paste(paste("/Users/grantsnitker/Desktop/WILDPLANTS/wild_plant_models/",paste(spNames[sp.i,1],spNames[sp.i,2],sep="_"),sep=""),"biomod_output/model.output",sep="/"))#setting the working directory to the models folder
myBiomodModelOut=get(load(modOut))#loading the model
setwd(paste(paste("/Users/grantsnitker/Desktop/WILDPLANTS/wild_plant_models/",paste(spNames[sp.i,1],spNames[sp.i,2],sep="_"),sep=""),"biomod_output",sep="/"))#as with the last two, all that really matters here is that the first string is correct

#this section of code will copy and paste the MaxEnt executable java file to a location where it can be used to reproect teh current species model
maxent.loc = getwd()
file.copy(from ="DIRECTORY TO MAXENT JAVA FILE", to = maxent.loc, overwrite = TRUE)
  
envStack=brick(MAP, MAT)#it is best to save your climate data as a brick (check out the R {raster} library), either as a .grd/gri, the native raster format for R, or as a .bil which is interchangeable with ArcGIS
envUnstack=unstack(envStack)#this and the next line are necessary to convert this to a raster stack instead of a brick
envStack=stack(envUnstack)
names(envStack) = c("MAP","MAT")
rm(envUnstack)#clearing memory
modProj=BIOMOD_Projection(modeling.output=myBiomodModelOut,new.env=envStack,proj.name='current',selected.models="all",binary.meth=NULL,compress="xz",clamping.mask=FALSE,keep.in.memory=F)#projecting the model (i.e. creating a map, takes awhile) and wrting it to file
projRast=(raster("model.output/proj_current/proj_current_model.output.grd")/1000)
reclass.projRast = reclassify(projRast, c(-Inf,threshold,0,threshold,Inf,1)) 
}
######## END WILD PLANT SPECIES RICHNESS FUNCTION  #######################


######## RUN FUNCTION  #######################
ptm <- proc.time()
wild_plant_result = brick(lapply(sp.number, function(i){wild.plants.function(i)}))
proc.time() - ptm

######## PROCESS AND EXPORT RESULTS  #######################
wild_plant_sum = calc(wild_plant_result,sum)
writeRaster(wild_plant_sum, filename = "OUTPUT DIRECTORY", format="GTiff", overwrite=TRUE)
