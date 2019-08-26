## LTVTP Anti-Correlation Script  v.16
## Developed by Grant Snitker; grant.snitker@uga.edu
## Dept of Crop and Soil Sciences
## University of Georgia
## Date: 08/26/2019

### Note: Replace data inputs (in CAPITAL LETTERS) with data for specified analytical region

######## BEGIN PREAMBLE ########
# Set the working directory.
setwd("PATH TO WORKING DIRECTORY")
# Load required packages
library("parallel")
library("plyr")
library("raster")
######## END PREAMBLE ########

######## BEGIN PARAMETERS ########
# Force Raster to load large rasters into memory
rasterOptions(chunksize=2e+07,maxmemory=2e+08)
# Import precipitation map for analytical region (the LTVTP project uses precip. data from PRISM Climate Group, Oregon State University, http://prism.oregonstate.edu )
prcp <- brick("RAW PRCP DATA") # rasterbrick of reconstructed precipitation values from the study region and temporal period of interest
prcp.mean = brick("RAW PRCP DATA MEAN") # rasterbrick of the mean reconstructed precipitation values from the study region and temporal period of interest
prcp.sd = brick("RAW PRCP DATA SD") # rasterbrick of the standard deviation of reconstructed precipitation values from the study region and temporal period of interest

### Create Site locations
Site1 = cellFromXY(prcp, c(0.000000, 0.000000)) # Enter the latitude and longitude of sites for the analysis here
Site2 = cellFromXY(prcp, c(0.000000, 0.000000))
Site3 = cellFromXY(prcp, c(0.000000, 0.000000))

# Create a list of all sites
all.sites = list(Site1, Site2, Site3)
# Output Type. Enter either "Matrix" or "Map" and reults will be exported in that format
Output = ""
######## END PARAMETERS ########

######## BEGIN ANTICOR FUNCTION ########
all.site.anticor = function(sites){
  
  # Parameters for run
  site = sites
  time.range = (1351:1450)# 1100:1250, 1251:1350, 1351:1450 this is the year range of interest from reconstruction
  
  # Total anti-corr function
  total.anticor = function(time.range){
    
    # Creating moving count for subsampling mean and sd
    count.sub <-count(time.range)
    count.sub = ((count.sub[1,1])-1)
    moving.years = time.range - count.sub
    
    # Create mean and sd for each cell
    prcp.sub <-subset(prcp,time.range)
    prcp.sub.mean <- subset(prcp.mean,moving.years)
    prcp.sub.sd <- subset(prcp.sd,moving.years)
    
    # Create upper and lower thresholds for each cell.  This will be used to assign anti-correlation
    upper.threshold = (prcp.sub.mean + (prcp.sub.sd * 0.5)) # upper threshold
    lower.threshold = (prcp.sub.mean - (prcp.sub.sd * 0.5)) #lower threshold
    
    # Positive Anti-correlation function
    positive.anticor = function(time.range){
      site.cell = extract(prcp.sub, site)
      focal.cell = site
      threshold =  extract(lower.threshold, focal.cell)
      output.map = prcp.sub
      output.map[] = NA
      
      compare.cells.pos <- function(prcp.map, output.map){
        if (site.cell < threshold) {output.map=overlay(prcp.map,upper.threshold, fun=function(x,y){return(x-y)})}
        else output.map[] = 0
        return(output.map)
      }
      out.pos = compare.cells.pos(prcp.sub,output.map)
      out.rc.pos <- reclassify(out.pos, c(-Inf,-0.001,0, 0,Inf,1))
      return(out.rc.pos)
    }
    
    # Negative Anti-correlation function
    negative.anticor = function(time.range){
      site.cell = extract(prcp.sub, site)
      focal.cell = site
      threshold =  extract(upper.threshold, focal.cell)
      output.map = prcp.sub
      output.map[] = NA
      
      compare.cells.neg <- function(prcp.map, output.map){
        if (site.cell > threshold) {output.map=overlay(prcp.map,lower.threshold, fun=function(x,y){return(x-y)})}
        else output.map[] = 0
        return(output.map)
      }
      out.neg = compare.cells.neg(prcp.sub,output.map)
      out.rc.neg <- reclassify(out.neg, c(-Inf,-0.001,1, 0,Inf,0))
      return(out.rc.neg) 
    }
    
    Pos.anticor = brick(lapply(time.range, function(i){positive.anticor(i)}))
    Pos.anticor.sum = (calc(Pos.anticor, sum))
    
    Neg.anticor = brick(lapply(time.range, function(i){negative.anticor(i)}))
    Neg.anticor.sum = (calc(Neg.anticor, sum))
    returnlist = list(Pos.anticor.sum, Neg.anticor.sum)
  }
  
  # Run total anti-corr function for year range and focal site
  total.anticor.result.nested = lapply(time.range, function(i){total.anticor(i)})
  
  
  # Process and un-nest results
  total.anticor.result = unlist(total.anticor.result.nested)
  count.total = count(time.range)
  Pos.anticor.sum <- calc(brick(total.anticor.result[1:length(total.anticor.result) %% 2 == 1]),sum)
  Neg.anticor.sum <- calc(brick(total.anticor.result[1:length(total.anticor.result) %% 2 == 0]),sum)
  Pos.anticor.per = round(Pos.anticor.sum/(sum(count.total[,2])),digits = 2)
  Neg.anticor.per = round(Neg.anticor.sum/(sum(count.total[,2])),digits = 2)
  Neg.anticor.per = round((Neg.anticor.sum/(sum(count.total[,2]))*1000),digits = -1)
  Combined_Anticor = (Pos.anticor.per + Neg.anticor.per)

  if(Output == "Map") {
    dwriteRaster(Combined_Anticor, filename=paste("OUTPUT PATH",paste(site,sep="_"), sep=""), format = "GTiff", overwrite=TRUE)}
  
  else {
  extractor = function(all.sites){
    ext = extract(Combined_Anticor, all.sites)
  }
  extracted.values = lapply(all.sites,extractor)
  }
}
######## END ANTICOR FUNCTION ########

######## RUN ANTICOR FUNCTION ########

### Create and Export Anticor results
ptm <- proc.time()
anticor.matrix.values = mclapply(sites, function(i){all.site.anticor(i)}, mc.cores = 4) # be sure to specific mc.cores to reflect the number of cores in the computer
proc.time() - ptm



##########################################################################################################
### If exporting results as a matrix, run the following lines of code to complete the processing and export
anticor.matrix.values.unlisted = unlist(anticor.matrix.values)
anticor.matrix = matrix(data = anticor.matrix.values.unlisted, nrow = 16, ncol = 16, byrow = T)# Adjust nrow and ncol to reflect the number of sites in the analysis.  For example, if there are 3 sites in the analysis, ncol and nrow should equal 3.
anticor.matrix.rast = raster(anticor.matrix)
writeRaster(anticor.matrix.rast, filename="OUTPUT PATH", format = "GTiff", overwrite=TRUE)
