
packrat::init("~/Dropbox/Gila Monster Data/GM_Study")

########################
## SPATIAL ANALYSES   ##
########################
# "./M69/Monsoon .csv"

# required packages
library(adehabitatHR) #for home range calculations
library(data.table) #manipulate S3 and S4 data tables
library(ggplot2) #for graphic output
library(ggfortify) #to allow ggplot2 to read spatial data
library(grid) #to add annotations to the output
# library(OpenStreetMap) #for obtaining raster images
library(pbapply) #needed for progress bar
library(plotly) #for interactive xy plot
library(rgdal) #for converting spatial data
library(sp) #for converting spatial data
library(rgeos)
# library(raster)
library(mapview)


########################
## PLOT ALL LIZARDS:
########################

#######################
# Plot spdf of all animal locations
GM.df<-read.csv("GM_Final_Data.csv")

# Using the function "SpatialPoints" we create an object of class SpatialPoints.
# We have to specify the coordinates, whereas the bbox is automatically generated.
GM.sp <- SpatialPoints(GM.df[,c("EASTING","NORTHING")], proj4string=CRS.SC)

# can easily plot it
plot(GM.sp, axes=TRUE)
# View(GM.df)
GM.spdf <- SpatialPointsDataFrame(coords=GM.df[,6:7], data=GM.df[,1:4], proj4string=CRS.SC)

plot(GM.spdf,axes=TRUE)

plot(GM.spdf, axes=TRUE, pch=21, cex=.2) # change point size

# 4.2 Advanced plotting with spplot()
spplot(GM.spdf, axes=TRUE, pch=21, cex=.2, zcol="YEAR", col.regions=rainbow(7))

# Advanced mapping with mapview
mapview(GM.spdf)

mapview(GM.spdf, cex=.2, zcol="LIZARDNUMBER", col.regions=rainbow(18), legend=F)

#########################################
## Plot all individuals for HR analyses

#load data with fake date and time
# data <- read.csv("GM_Final_Data.csv")
# View(data)
# 
# #interactive examination of data points for outliers with Plotly
# p <- ggplot() + geom_point(data=data, aes(EASTING,NORTHING, color=LIZARDNUMBER)) +
#   labs(x="Easting", y="Northing")
# ggplotly(p)
# 
# #split into multiple files for individuals
# lapply(split(data, data$LIZARDNUMBER), 
#        function(x)write.csv(x, file = paste(x$LIZARDNUMBER[1],".csv"), row.names = FALSE))
# 
# #create list of individual files created in previous step #edit pattern based on individual ids
# files <- list.files(path = ".", pattern = "[MF]+[0-9]", full.names = TRUE)
# 
# # Set proj.string
CRS.SC<-CRS("+proj=utm +zone=12 +ellps=WGS84 +units=m +no_defs")

# #creating spatial data frame for all points
# x <- as.data.frame(data$EASTING)
# y <- as.data.frame(data$NORTHING)
# xy <- c(x,y)
# data.proj <- SpatialPointsDataFrame(xy,data, proj4string = CRS.SC)
# 
# #creating homerange in adehabitat for all points
# xy <- SpatialPoints(data.proj@coords)
# mcp.out <- mcp(xy, percent=100, unout="ha")
# #plot(xy)
# #plot(mcp.out)
# mcp_area <- as.data.frame(mcp.out@data$area)
# colnames(mcp_area) <- "Hectares"
# write.table(mcp_area, "MCP_Area.csv", sep = ",", row.names = TRUE)
# 
# #KDE creation for all points
# kde<-kernelUD(xy, h="href", kern="bivnorm", grid=1000)
# ver <- getverticeshr(kde, 95)
# ver$area
# kde@h$h
# #plot(data.proj@coords)
# #plot(ver)
# kde_area <- as.data.frame(ver$area)
# colnames(kde_area) <- "Hectares"
# write.table(kde_area, "KDE_Area.csv", sep = ",", row.names = TRUE)
# 
# #mcp plot for all points
# mcp.points <- cbind((data.frame(xy)),data$LIZARDNUMBER)
# colnames(mcp.points) <- c("x","y", "lizardnumber")
# mcp.poly <- fortify(mcp.out, region = "id")
# mcp.plot <- ggplot()+
#   geom_polygon(data=mcp.poly, aes(x=mcp.poly$long, y=mcp.poly$lat))+
#   geom_point(data=mcp.points, aes(x=x, y=y,color = lizardnumber))
# mcp.plot
# 
# #kde plot for all points
# kde.points <- cbind((data.frame(data.proj@coords)),data$LIZARDNUMBER)
# colnames(kde.points) <- c("x","y","lizardnumber")
# kde.poly <- fortify(ver, region = "id")
# kde.plot <- ggplot()+
#   geom_polygon(data=kde.poly, aes(x=kde.poly$long, y=kde.poly$lat))+
#   geom_point(data=kde.points, aes(x=x, y=y, color = lizardnumber))
# kde.plot

########################## INDIVIDAL ANALYSES #########################
# BY YEAR
# MCP without raster
# Code set to run individual lizard by ID. Files for by year analyses are set in ID 
# subfolders.
#####################

# BREAK DOWN INDIVIUAL ANIMAL DATA BY YEAR and or SEASON, *NOTE MUST GO INTO DROPBOX FOR EACH ONE AND 
# MOVE "YEAR" FILES TO THEIR INDIVIDUAL FOLDERS
# setwd("~/Dropbox/Gila Monster Data/GM_Study/F104")
# library(readr)
# M112<-read_csv("M112 .csv")
# View(F66)
# 
# # split into multiple files for individual for YEAR
# lapply(split(M112, M112$YEAR),
#        function(x)write.csv(x, file = paste(x$YEAR[3],".csv"), row.names = FALSE))
# 
# # split into multiple files for individual for SEASON
# lapply(split(F104, F104$SEASON),
#        function(x)write.csv(x, file = paste(x$SEASON[4],".csv"), row.names = FALSE))

# remove unwanted objects
# rm(F36,M215) 

###################################################
##        INDIVIDUAL SPATIAL ANALYSIS            ## 
###################################################

###############################
# Plot spdf of animal locations

M69.df<-read.csv("./M69/M69 .csv")
M67.df<-read.csv("./M67/M67 .csv")
M255.df<-read.csv("./M255/M255 .csv")
M215.df<-read.csv("./M215/M215 .csv")
M14.df<-read.csv("./M14/M14 .csv")
M119.df<-read.csv("./M119/M119 .csv")
M112.df<-read.csv("./M112/M112 .csv")
View(M215.df)
# Using the function "SpatialPoints" we create an object of class SpatialPoints.
# We have to specify the coordinates, whereas the bbox is automatically generated.
# M67.sp <- SpatialPoints(M67.df[,c("EASTING","NORTHING")], proj4string=CRS.SC)

# can easily plot it
# plot(M67.sp, axes=TRUE)
# View(GM.df)
M69 <- SpatialPointsDataFrame(coords=M69.df[,6:7], data=M69.df[,1:4], proj4string=CRS.SC)
M67 <- SpatialPointsDataFrame(coords=M67.df[,6:7], data=M67.df[,1:4], proj4string=CRS.SC)
M255 <- SpatialPointsDataFrame(coords=M255.df[,7:8], data=M255.df[,1:5], proj4string=CRS.SC)
M215 <- SpatialPointsDataFrame(coords=M215.df[,7:8], data=M215.df[,1:5], proj4string=CRS.SC)
M14 <- SpatialPointsDataFrame(coords=M14.df[,7:8], data=M14.df[,1:5], proj4string=CRS.SC)
M119 <- SpatialPointsDataFrame(coords=M119.df[,7:8], data=M119.df[,1:5], proj4string=CRS.SC)
M112 <- SpatialPointsDataFrame(coords=M112.df[,7:8], data=M112.df[,1:5], proj4string=CRS.SC)

# plot(M67,axes=TRUE)

plot(M69, axes=TRUE, pch=21, cex=.5) # change point size

# 4.2 Advanced plotting with spplot()
spplot(M69, axes=TRUE, pch=21, cex=.5, zcol="YEAR", col.regions=rainbow(3))

# Advanced mapping with mapview
mapview(M67)

# mapview(M67, cex=.5, zcol="LIZARDNUMBER", col.regions=rainbow(1), legend=F)


#########################
##        MCP          ## 
#########################

# Set projection forestring
#CRS.SC<-CRS("+proj=utm +zone=12 +ellps=WGS84 +units=m +no_defs")

# SET WORKING DIRECTORY TO INDIDUAL SUBFOLDER
library(readr)
# setwd("~/Dropbox/Gila Monster Data/GM_Study/M67")
# F104<-read_csv("F104 .csv")
# View(F104)

# Function for running animal MCP:
MCP<-mcp_analysis("./F114/2007 .csv", percentage=95)
MCP<-mcp_analysis("./F36/2010 .csv", percentage=95)
summary(MCP)

###########################################################################################
## looping function for mcp ******BY YEAR******

mcp_analysis <- function(filename, percentage){
  data <- read.csv(file = filename)
  x <- as.data.frame(data$EASTING)
  y <- as.data.frame(data$NORTHING)
  xy <- c(x,y)
  data.proj <- SpatialPointsDataFrame(xy,data, proj4string = CRS.SC)
  xy <- SpatialPoints(data.proj@coords)
  mcp.out <- mcp(xy, percentage, unout="ha")
  area <- as.data.frame(round(mcp.out@data$area,4))
  .rowNamesDF(area, make.names=TRUE) <- data$YEAR
  write.table(area,file="MCP_Hectares.csv",
              append=TRUE,sep=",", col.names=FALSE, row.names=TRUE)
  mcp.points <- cbind((data.frame(xy)),data$YEAR)
  colnames(mcp.points) <- c("x","y", "year")
  mcp.poly <- fortify(mcp.out, region = "id")
  units <- grid.text(paste(round(mcp.out@data$area,2)," ha"), x=0.9,  y=0.95,
                     gp=gpar(fontface=4, cex=0.9), draw = FALSE)
  mcp.plot <- ggplot() +
    geom_polygon(data=mcp.poly, aes(x=mcp.poly$long, y=mcp.poly$lat), alpha=0.5) +
    geom_point(data=mcp.points, aes(x=x, y=y)) + theme_bw() +
    labs(x="Easting (m)", y="Northing (m)", title=mcp.points$year) +
    theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5)) +
    annotation_custom(units)
  mcp.plot
}

## looping function for mcp ******BY YEAR******
###########################################################################################

#individual mcp run
# mcp_analysis("./2008 .csv", percentage=50)
# #run all individuals for mcp, takes 12sec
# lapply(files,mcp_analysis)
# #pblapply(files, mcp_analysis) #runs with progressbar




####################################################
### *****WORKING***** MCP Polygon *by year  
####################################################

M67_MCP<-mcp_analysis.POLY('./M67/M67 .csv', percentage= 100)
M69_MCP<-mcp_analysis.POLY('./M67/M67 .csv', percentage= 100)
M255_MCP<-mcp_analysis.POLY('./M67/M67 .csv', percentage= 100)
M215_MCP<-mcp_analysis.POLY('./M67/M67 .csv', percentage= 100)
M14_MCP<-mcp_analysis.POLY('./M67/M67 .csv', percentage= 100)
M119_MCP<-mcp_analysis.POLY('./M67/M67 .csv', percentage= 100)
M112_MCP<-mcp_analysis.POLY('./M67/M67 .csv', percentage= 100)

M67_MCP

######################################################################
##  ***USE THIS*** FOR MCP POLYGON OUTPUTS FOR INTERSECT CALCULATIONS:


## File name = X_MCP

mapviewOptions(basemaps = c("OpenStreetMap","Esri.WorldImagery","OpenTopoMap"),
               na.color = "magenta",
               layers.control.pos = "topleft")

M112_MCP<-mcp_analysis.POLY('./M112/M112 .csv', percentage= 100)
mapView(M112_MCP, zcol="id")

############################################################################################
#################### function #################### function ################################

# mcp_analysis.POLY <- function(filename, percentage){
#   data <- read.csv(file = filename,stringsAsFactors = FALSE)
#   data.sp <- data[, c("LIZARDNUMBER", "EASTING", "NORTHING")]
#   coordinates(data.sp) <- c("EASTING", "NORTHING")
#   proj4string(data.sp) <- CRS.SC
#   mcp_out <- mcp(data.sp, percentage, unout="ha")
# }

#################### function  ################### function ################################
############################################################################################





########################################################
##                    KDE ANALYSES                    ##
##                  *** WORKING ***                   ##
########################################################


# Set projection forestring, already done previously, but can change if needed:
# CRS.SC<-CRS("+proj=utm +zone=12 +ellps=WGS84 +units=m +no_defs") 

# Function for running animal KDE:
kde_analysis.href.plot("./F66/2011 .csv", percentage=95)
#lapply(files,kde_analysis.href.plot, percentage=95)


#########################################
##  looping function for kde ****href****

################# looping function################### looping function ################

kde_analysis.href.plot <- function(filename, percentage){
  data <- read.csv(file = filename)
  x <- as.data.frame(data$EASTING)
  y <- as.data.frame(data$NORTHING)
  xy <- c(x,y)
  data.proj <- SpatialPointsDataFrame(xy,data, proj4string = CRS.SC)
  xy <- SpatialPoints(data.proj@coords)
  kde<-kernelUD(xy, h="href", kern="bivnorm", grid=1000)
  ver <- getverticeshr(kde, percentage)
  area <- as.data.frame(round(ver$area,4))
  .rowNamesDF(area, make.names=TRUE) <- data$LIZARDNUMBER
  write.table(area,file="KDE_Hectares.csv",
              append=TRUE,sep=",", col.names=FALSE, row.names=TRUE)
  kde.points <- cbind((data.frame(data.proj@coords)),data$LIZARDNUMBER)
  colnames(kde.points) <- c("x","y","lizardnumber")
  kde.poly <- fortify(ver, region = "id")
  units <- grid.text(paste(round(ver$area,2)," ha"), x=0.9,  y=0.95,
                     gp=gpar(fontface=4, cex=0.9), draw = FALSE)
  kde.plot <- ggplot() +
    geom_polygon(data=kde.poly, aes(x=kde.poly$long, y=kde.poly$lat), alpha = 0.5) +
    geom_point(data=kde.points, aes(x=x, y=y)) + theme_bw() +
    labs(x="Easting (m)", y="Northing (m)", title=kde.points$lizardnumber) +
    theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5)) +
    annotation_custom(units)
  kde.plot
}

# individual kde with ***href*** run Plot, with change to "percentage"
#kde_analysis.href.plot("./2008 .csv", percentage=95)
# run all individuals for kde, takes 3min
#lapply(files,kde_analysis, percentage=95)     # **NOT TESTED**
#pblapply(files, kde_analysis, percentage=95)  #runs with progressbar  **NOT TESTED**


######################################################################
######  kde POLYGON  href, by year or by all years
######################################################################

# set working directory for individual animals if needed:
#setwd("~/Dropbox/Gila Monster Data/GM_Study/M69")

M112_KDEPoly<-kde_analysis.href.polygon("./M112 .csv", percentage= 95)
mapView(M112_KDEPoly)

F36_95<-kde_analysis.href.polygon("./F36 .csv", percentage= 95)
#plot(F36_95)
mapView(F36_95)
F36_95
## change polygon opacity using alpha.regions=x in mapview:
mapView(F36_95, alpha.regions=0.3)


###############################################
##  Mapview with multiple polygons ***WORKING***:

mapView(M69_95+F147_95)

library(plainview)
# mapview w list of objects
#mapview(M69_95,alpha.regions=0.3) + mapview(F147_95,alpha.regions=0.2)

# Have to specify each attribute needed then reassign to an abject using mapview()
m1 = mapview(M69_95, zcol = "id", col.regions = c("blue"),alpha.regions=0.3)
m2 = mapview(F147_95, zcol = "id", col.regions = c("red"),alpha.regions=0.3)
m1
m2
m3 = mapview(M119_95, zcol = "id", col.regions = c("blue"),alpha.regions=0.3)
m4 = mapview(F36_95, zcol = "id", col.regions = c("red"),alpha.regions=0.3)

# add each newly assign polygon to eachother for n layers
m1+m2
m3+m4

############################################################################################
####################  Looping function #################### Looping function ###############

# kde_analysis.href.polygon <- function(filename, percentage){
#   data <- read.csv(file = filename)
#   x <- as.data.frame(data$EASTING)
#   y <- as.data.frame(data$NORTHING)
#   xy <- c(x,y)
#   data.proj <- SpatialPointsDataFrame(xy,data, proj4string = CRS.SC)
#   xy <- SpatialPoints(data.proj@coords)
#   kde<-kernelUD(xy, h="href", kern="bivnorm", grid=1000)
#   ver <- getverticeshr(kde, percentage)
#   ver@proj4string<-CRS.SC
#   area <- as.data.frame(round(ver$area,4))
#   .rowNamesDF(area, make.names=TRUE) <- data$YEAR
#   write.table(area,file="KDE_Hectares.csv",
#               append=TRUE,sep=",", col.names=FALSE, row.names=TRUE)
#   kde.points <- cbind((data.frame(data.proj@coords)),data$YEAR)
#   colnames(kde.points) <- c("x","y","year")
#   kde.poly <- fortify(ver, region = "id")
#   units <- grid.text(paste(round(ver$area,2)," ha"), x=0.9,  y=0.95,
#                      gp=gpar(fontface=4, cex=0.9), draw = FALSE)
#   ver
# }

####################  Looping function  ################### Looping function ###############
############################################################################################

# create the polygon shapefile and plot:
#setwd("~/Dropbox/Gila Monster Data/GM_Study/F36")

# library(mapview)
# F36_95<-kde_analysis.href.polygon("./F36 .csv", percentage= 95)
# #plot(F36_95)
# mapview(F36_95)

# # For all Animals
# ALL_95<-kde_analysis.href.polygon("./SC_Year_Dis.csv", percentage= 95)
# #plot(ALL_95)
# mapview(ALL_95)

# remove unwanted objects
#rm(F147_95,F2007,M119_95,M69_95)

############################
###### Create raster of UD

#setwd("~/Dropbox/Gila Monster Data/GM_Study/F36")

# kde_analysis.href.raster <- function(filename){
#   data <- read.csv(file = filename)
#   x <- as.data.frame(data$EASTING)
#   y <- as.data.frame(data$NORTHING)
#   xy <- c(x,y)
#   data.proj <- SpatialPointsDataFrame(xy,data, proj4string = CRS.SC)
#   xy <- SpatialPoints(data.proj@coords)
#   kde<-kernelUD(xy, h="href", kern="bivnorm", grid=1000)
#   kde<-as(kde, "SpatialGridDataFrame")
#   kde@proj4string<- CRS.SC
#   kde
# }

# Create raster object, can plot it and/or overlay it over imagery in mapview
library(mapview)
F36.raster<-kde_analysis.href.raster("./2008 .csv")
mapview(F36.raster)
# plot(F36.raster)

# REMOVE UNWANTED FILES
#rm(M69_POLY,M69.95)

#################################################################
##                Temporal Trajectory/Distance                 ##
#################################################################

#*********WORKING**********

# #trajectory analysis and distance over time

## Individual animal should already be assighned to an object from above. This set should
## also be across all years.  For analysis of each year, you will just use the "year .csv"
## file within the subfolder.

# View(F36)

# trajectory analysis:

#pblapply(files, traj_analysis)
traj_analysis("M67 .csv")
traj_analysis("2007 .csv")

############################################################################################
####################  Looping function #################### Looping function ###############

# traj_analysis <- function(filename){
#   relocs_data <- read.csv(file = filename)
#   relocs <- as.ltraj(cbind(relocs_data$EASTING, relocs_data$NORTHING),id=relocs_data$LIZARDNUMBER, typeII = FALSE, date=NULL)
#   relocs.df <- ld(relocs)
#   relocs_dist <- as.data.frame(sum(sapply(relocs.df$dist, sum, na.rm=TRUE)))
#   colnames(relocs_dist) <- "Total Distance"
#   name <- relocs.df$id[1]
#   row.names(relocs_dist) <- name
#   relocs_units <- grid.text(paste(round(relocs_dist,2),"m"), x=0.9, y=0.9,
#                             gp=gpar(fontface=3, col="black", cex=0.9), draw = FALSE)
#   reloc.plot <- ggplot() + theme_classic() + geom_path(data=relocs.df, aes(x=x,y=y), linetype = "dashed", colour = "red",
#                                                        arrow = arrow(length=unit(.5,"cm"), angle = 20, ends="last", type = "closed")) +
#     geom_point(data=relocs.df, aes(x=x, y=y)) + geom_point(data=relocs.df, aes(x=x[1],
#                                                                                y=y[1]), size = 3, color = "darkgreen", pch=0) +
#     labs(x="Easting (m)", y="Northing (m)", title=relocs.df$id[1]) +
#     theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5)) +
#     annotation_custom(relocs_units)
#   reloc.plot
# }

####################  Looping function  ################### Looping function ###############
############################################################################################



##########################################
# TRAJECTORY ANALYSIS, VECTOR FOR MAPVIEW:

# ************ WOORKING ******************** WORKING ************************* WORKING *****

M67_Lines <- traj_analysis.lines("M67 .csv")
mapView(M67_Lines)


############################################################################################
####################  Looping function #################### Looping function ###############

traj_analysis.lines <- function(filename){
  relocs_data <- read.csv(file = filename)
  relocs <- as.ltraj(cbind(relocs_data$EASTING, relocs_data$NORTHING),id=relocs_data$LIZARDNUMBER, 
                     typeII = FALSE, date=NULL, proj4string = CRS.SC)
  relocs.df <- ld(relocs)
  relocs_dist <- as.data.frame(sum(sapply(relocs.df$dist, sum, na.rm=TRUE)))
  colnames(relocs_dist) <- "Total Distance"
  name <- relocs.df$id[1]
  row.names(relocs_dist) <- name
  relocs_units <- grid.text(paste(round(relocs_dist,2),"m"), x=0.9, y=0.9,
                            gp=gpar(fontface=3, col="black", cex=0.9), draw = FALSE)
  relocs
}

?ld()
####################  Looping function  ################### Looping function ###############
############################################################################################


##############################
# distance over time analysis:

#pblapply(files, dist_analysis)
dist_analysis("F36 .csv")
dist_analysis("2007 .csv")

#remove individual files after analysis, this works so be careful!
#unlink(files, recursive = FALSE)

# dist_analysis <- function(filename){
#   relocs_data <- read.csv(file = filename)
#   relocs <- as.ltraj(cbind(relocs_data$EASTING, relocs_data$NORTHING),id=relocs_data$LIZARDNUMBER, typeII = FALSE, date=NULL)
#   relocs.df <- ld(relocs)
#   relocs_dist <- as.data.frame(sum(sapply(relocs.df$dist, sum, na.rm=TRUE)))
#   colnames(relocs_dist) <- "Total Distance"
#   name <- relocs.df$id[1]
#   row.names(relocs_dist) <- name
#   write.table(relocs_dist,file="reloc_dist.csv",
#               append=TRUE,sep=",", col.names=FALSE, row.names=TRUE)
#   dist.plot
# }



############################################
##          TRACKING INTENSITY            ##
############################################


library(tidyverse)
library(lubridate)
library(plotly)
library(scales)


TrackData <- read.csv("./SC_movement.csv")
View(TrackData)

TrackData$DATE <- mdy(TrackData$DATE) #use lubridate to specify incoming date format; 
#tell lubridate that the current format is "mdy", it then converts column to "ymd" format
str(TrackData) #look at DATE format

pT <- ggplot(TrackData, aes(DATE, LIZARDNUMBER)) +
  geom_count (color = "blue") + #set symbol size proportional to # of overlapping observations (same day)
  scale_x_date(date_breaks = "12 month", labels = date_format("%Y"))
pT

#####################################
##          Bootstrap MCP          ##
#####################################

getClass("Spatial")
getClass("SpatialPoints")

# An object of class SpatialPoints is a Spatial with one more slot:
# 3) coords: a matrix containing the coordinates of the points

# upload a dataframe with XY coordinates

setwd("~/Dropbox/Gila Monster Data/GM_Study/M255") # please set your working directory
M255.df<-read.csv("M255 .csv")

# # have a look at it
# head(F36.df) # First few rows
# summary(F36.df) # Statistical summaries
# str(F36.df) # Aspects of its structure


# Using the function "SpatialPoints" we create an object of class SpatialPoints.
# We have to specify the coordinates, whereas the bbox is automatically generated.
M255.sp <- SpatialPoints(M255.df[,c("EASTING","NORTHING")], proj4string=CRS.SC)
# summary(F36.sp)
# bbox(F36.sp)
# proj4string(F66.sp)
coordinates(M255.sp) # to have a look at the coordinates slot

plot(M255.sp, axes=TRUE)

library(move)
# ?move

hrBootstrap(x=M255.sp, rep=100, unin='m', unout='ha')

# rm(M255.df)

################################
##      SPATIAL ANALYSES      ##
##      MCP MAPS/OVERLAP      ##
################################

library(rgeos)

# creating interactive map of male vs. female MCPs:
# mapview(M67_MCP,legend=F)+mapview(M69_MCP,legend=F)+mapview(M255_MCP,legend=F)+
#   mapview(M215_MCP,legend=F)+mapview(M14_MCP,legend=F)+mapview(M119_MCP,legend=F)+
#   mapview(M112_MCP,legend=F)+mapview(F66_MCP,legend=F, zcol = "id", col.regions = c("red"),alpha.regions=0.3)+
#   mapview(F36_MCP,legend=F, zcol = "id", col.regions = c("red"),alpha.regions=0.3)+
#   mapview(F252_MCP,legend=F, zcol = "id", col.regions = c("red"),alpha.regions=0.3)+ 
#   mapview(F214_MCP,legend=F, zcol = "id", col.regions = c("red"),alpha.regions=0.3)+
#   mapview(F200_MCP,legend=F, zcol = "id", col.regions = c("red"),alpha.regions=0.3)+
#   mapview(F147_MCP,legend=F, zcol = "id", col.regions = c("red"),alpha.regions=0.3)+
#   mapview(F146_MCP,legend=F, zcol = "id", col.regions = c("red"),alpha.regions=0.3)+
#   mapview(F137_MCP,legend=F, zcol = "id", col.regions = c("red"),alpha.regions=0.3)+
#   mapview(F135_MCP,legend=F, zcol = "id", col.regions = c("red"),alpha.regions=0.3)+
#   mapview(F114_MCP,legend=F, zcol = "id", col.regions = c("red"),alpha.regions=0.3)+
#   mapview(F104_MCP,legend=F, zcol = "id", col.regions = c("red"),alpha.regions=0.3)

## Shortcut is to rbind() males together and females together then use mapview to map them:

Male.MCP <- rbind(M67_MCP,M69_MCP,M255_MCP,M215_MCP,M14_MCP,M119_MCP,M112_MCP)
mapview(Male.MCP)
Female.MCP <- rbind(F66_MCP,F36_MCP,F252_MCP,F214_MCP,F200_MCP,F147_MCP,F146_MCP,F137_MCP,
                    F135_MCP,F114_MCP,F104_MCP)
# mapview(Female.MCP)

mapview(Male.MCP, legend=F, zcol="id", col.regions = c("blue"), alpha.regions=0.3) + 
  mapview(Female.MCP, legend=F, zcol = "id", col.regions = c("red"), alpha.regions=0.3)

##############################
## MCP OVERLAP BETWEEN SEXES:

Male.MCP
Female.MCP

MCP_Intersect<-gIntersection(Male.MCP, Female.MCP,byid=T)
MCP_Intersect$area<-gArea(MCP_Intersect, byid=T)/10000
MCP_Intersect
mapView(MCP_Intersect,legend=F, col.regions = c("purple"), alpha.regions=0.3)

## A table:
kable(MCP_Intersect, format = "pandoc", caption = 'Home Range Overlap by Sex')

#################################
## NET MCP OVERLAP WITHIN SEXES :

F_OL1<-rbind(F36_MCP,F146_MCP)
F_OL2<-rbind(F66_MCP,F146_MCP)
F_OL3<-rbind(F104_MCP,F137_MCP,F147_MCP)
F_OL4<-rbind(F36_MCP,F66_MCP)
F_OL5<-rbind(F137_MCP,F135_MCP)

Female_Intersect_1<-gIntersection(F66_MCP,F_OL1,byid=F)
Female_Intersect_1$area<-gArea(Female_Intersect_1, byid=T)/10000
Female_Intersect_1
mapView(Female_Intersect_1,legend=F, col.regions = c("red"), alpha.regions=0.3)

Female_Intersect_2<-gIntersection(F36_MCP,F_OL2,byid=F)
Female_Intersect_2$area<-gArea(Female_Intersect_2, byid=T)/10000
Female_Intersect_2
mapView(Female_Intersect_2,legend=F, col.regions = c("red"), alpha.regions=0.3)

Female_Intersect_2a<-gIntersection(F104_MCP,F135_MCP,byid=F)
Female_Intersect_2a$area<-gArea(Female_Intersect_2a, byid=T)/10000
Female_Intersect_2a
mapView(Female_Intersect_2a,legend=F, col.regions = c("red"), alpha.regions=0.3)

Female_Intersect_3<-gIntersection(F135_MCP,F_OL3,byid=F)
Female_Intersect_3$area<-gArea(Female_Intersect_3, byid=T)/10000
Female_Intersect_3
mapView(Female_Intersect_3,legend=F, col.regions = c("red"), alpha.regions=0.3)

Female_Intersect_4<-gIntersection(F146_MCP,F_OL4,byid=F)
Female_Intersect_4$area<-gArea(Female_Intersect_4, byid=T)/10000
Female_Intersect_4
mapView(Female_Intersect_4,legend=F, col.regions = c("red"), alpha.regions=0.3)

Female_Intersect_5<-gIntersection(F147_MCP,F_OL5,byid=F)
Female_Intersect_5$area<-gArea(Female_Intersect_5, byid=T)/10000
Female_Intersect_5
mapView(Female_Intersect_5,legend=F, col.regions = c("red"), alpha.regions=0.3)

## M215:M119
Male_Intersect_1<-gIntersection(M215_MCP, M119_MCP,byid=T)
Male_Intersect_1$area<-gArea(Male_Intersect_1, byid=T)/10000
Male_Intersect_1
mapView(Male_Intersect_1,legend=F, col.regions = c("blue"), alpha.regions=0.3)


## M14:M69
Male_Intersect_2<-gIntersection(M14_MCP,byid=T)
Male_Intersect_2$area<-gArea(Male_Intersect_2, byid=T)/10000
Male_Intersect_2
mapView(Male_Intersect_2,legend=F, col.regions = c("blue"), alpha.regions=0.3)

Female_Intersect_1<-gIntersection(F36_MCP, F146_MCP, byid=T)
Female_Intersect_1$area<-gArea(Female_Intersect_1, byid=T)/10000
Female_Intersect_1
mapView(Female_Intersect_1,legend=F, col.regions = c("blue"), alpha.regions=0.3)

#############################
## NET OVERLAP BETWEEN SEXES:

MF_OL1<-rbind(F135_MCP,F137_MCP,F147_MCP,F146_MCP,F66_MCP)

Net_Inter_1<-gIntersection(M69_MCP,MF_OL1,byid=F)
Net_Inter_1$area<-gArea(Net_Inter_1, byid=T)/10000
Net_Inter_1
mapView(Net_Inter_1,legend=F, col.regions = c("green"), alpha.regions=0.3)

## TOTAL NET OVERLAP OF MALES:FEMALES:
OL_Complex1<-rbind(F66_MCP,F146_MCP,M119_MCP,M215_MCP)
OL_Complex2<-rbind(F66_MCP,F146_MCP,M14_MCP,F147_MCP,F137_MCP,F135_MCP)
OL_Complex3<-rbind(M69_MCP,F147_MCP,F137_MCP,F104_MCP)
OL_Complex4<-rbind(F146_MCP,M69_MCP,F147_MCP)
OL_Complex5<-rbind(F137_MCP,F135_MCP,M69_MCP,M14_MCP,M67_MCP)
OL_Complex6<-rbind(F36_MCP,F66_MCP,M69_MCP,M14_MCP)
OL_Complex7<-rbind(F36_MCP,F146_MCP,M69_MCP)
OL_Complex8<-rbind(M215_MCP,F36_MCP)
OL_Complex9<-rbind(M119_MCP,F36_MCP)

Net_Intersect_1<-gIntersection(F36_MCP,OL_Complex1,byid=F)
Net_Intersect_1$area<-gArea(Net_Intersect_1, byid=T)/10000
Net_Intersect_1
mapView(Net_Intersect_1,legend=F, col.regions = c("green"), alpha.regions=0.3)

Net_Intersect_2<-gIntersection(M69_MCP,OL_Complex2,byid=F)
Net_Intersect_2$area<-gArea(Net_Intersect_2, byid=T)/10000
Net_Intersect_2
mapView(Net_Intersect_2,legend=F, col.regions = c("green"), alpha.regions=0.3)

Net_Intersect_3<-gIntersection(F135_MCP,OL_Complex3,byid=F)
Net_Intersect_3$area<-gArea(Net_Intersect_3, byid=T)/10000
Net_Intersect_3
mapView(Net_Intersect_3,legend=F, col.regions = c("green"), alpha.regions=0.3)

Net_Intersect_4<-gIntersection(M14_MCP,OL_Complex4,byid=F)
Net_Intersect_4$area<-gArea(Net_Intersect_4, byid=T)/10000
Net_Intersect_4
mapView(Net_Intersect_4,legend=F, col.regions = c("green"), alpha.regions=0.3)

Net_Intersect_5<-gIntersection(F147_MCP,OL_Complex5,byid=F)
Net_Intersect_5$area<-gArea(Net_Intersect_5, byid=T)/10000
Net_Intersect_5
mapView(Net_Intersect_5,legend=F, col.regions = c("green"), alpha.regions=0.3)

Net_Intersect_6<-gIntersection(F146_MCP,OL_Complex6,byid=F)
Net_Intersect_6$area<-gArea(Net_Intersect_6, byid=T)/10000
Net_Intersect_6
mapView(Net_Intersect_6,legend=F, col.regions = c("green"), alpha.regions=0.3)

Net_Intersect_7<-gIntersection(F66_MCP,OL_Complex7,byid=F)
Net_Intersect_7$area<-gArea(Net_Intersect_7, byid=T)/10000
Net_Intersect_7
mapView(Net_Intersect_7,legend=F, col.regions = c("green"), alpha.regions=0.3)

Net_Intersect_8<-gIntersection(M119_MCP,OL_Complex8,byid=F)
Net_Intersect_8$area<-gArea(Net_Intersect_8, byid=T)/10000
Net_Intersect_8
mapView(Net_Intersect_8,legend=F, col.regions = c("green"), alpha.regions=0.3)

Net_Intersect_9<-gIntersection(M215_MCP,OL_Complex9,byid=F)
Net_Intersect_9$area<-gArea(Net_Intersect_9, byid=T)/10000
Net_Intersect_9
mapView(Net_Intersect_9,legend=F, col.regions = c("green"), alpha.regions=0.3)



####################################
##        SPATIAL ANALYSES        ##
##     YEARLY/SEASONAL SHIFTS     ##
####################################

###########################################################################################
##                                       YEAR                                            ##


##########################################
## MCP Polygons for analyzing HR shifts         

# Plot spdf of animal locations
F66.df<-read.csv("./F66/F66 .csv")

# Using the function "SpatialPoints" we create an object of class SpatialPoints.
# We have to specify the coordinates, whereas the bbox is automatically generated.
F66.sp <- SpatialPoints(F66.df[,c("EASTING","NORTHING")], proj4string=CRS.SC)
# rm(F104.sp)

# can easily plot it
# plot(M67.sp, axes=TRUE)
# View(GM.df)
F66.spdf <- SpatialPointsDataFrame(coords=F66.df[,7:8], data=F66.df[,1:5], 
                                    proj4string=CRS.SC)
## Plot spdf:
# plot(F104.spdf, axes=TRUE, pch=21, cex=.5) # change point size

# Advanced plotting with spplot()
# spplot(F104.spdf, axes=TRUE, pch=21, cex=.5, zcol="SEASON", col.regions=rainbow(5))

## Create MCP polygons by year:
F66_mcp.08<-mcp_analysis.POLY("./2008 .csv", percentage= 100)

## Plot polygons:
plot(F104_mcp.08+F104_mcp.09,axes=TRUE)
plot(F104_mcp.08+F104_mcp.09,axes=TRUE,title(main = "F104 Home Range Shift"))

## Plot Points
# spplot(F104.spdf, axes=TRUE, pch=21, cex=.5, zcol="SEASON", col.regions=rainbow(5))

############################################################
## layered map with spplot() ******* WORKING *******

## bind polygons together
F104_mcp.shift<-rbind(F104_mcp.08,F104_mcp.09)
F66_mcp.shift<-rbind(F66_mcp.08,F66_mcp.09,F66_mcp.10)


# Create a layer-list
points.layer <- list("sp.points", F104.spdf,axes=TRUE, pch=21, cex=.5, col = "black")
polygons.layer <- list("sp.polygons", F104_mcp.shift, col = "blue", lwd=2)

points.layer.F66 <- list("sp.points", F66.spdf,axes=TRUE,  zcol="SEASON",pch=21, cex=.5, col = "black")
polygons.layer.F66 <- list("sp.polygons", F66_mcp.shift, col = "black", lwd=2)
# F66.SpPoly.layer <- rbind(points.layer.F66,polygons.layer.F66)

## Plot with layer-list passed to sp.layout`:
spplot(F104_mcp.shift, zcol="id", sp.layout = points.layer)
spplot(F104_mcp.shift, col="red", zcol="id")

spplot(F104.spdf, zcol="SEASON", sp.layout = polygons.layer)
spplot(F66.spdf, zcol="SEASON", sp.layout = polygons.layer.F66)



# You can also use the sp.layout option to add other things, like compass arrows or scales.
# You have to set the "offset", which defines the location of the 
# bottom left hand corner of the object you are working on.

F104_mcp.shift@bbox
# let put these in the right position, and add the text for scale
scale <- list("SpatialPolygonsRescale", layout.scale.bar(),fill = c("transparent","black"),
              offset = c(501330, 3591975), scale=100)
text1 <- list("sp.text", c(501330, 3591990), "0")
text2 <- list("sp.text", c(501430,  3591990), "100 m")

spplot(F104.spdf,zcol="SEASON", sp.layout = list(polygons.layer,scale,text1,text2))

########################### WORKING GGPLOT OF MULTI POLYGONS ###############################

## Create MCP polygons by YEAR:
M215_mcp.11<-mcp_analysis.POLY("./M215/2011 .csv", percentage= 100)
M215_mcp.12<-mcp_analysis.POLY("./M215/2012 .csv", percentage= 100)
F104_mcp.08<-mcp_analysis.POLY("./F104/2008 .csv", percentage= 100)
F104_mcp.09<-mcp_analysis.POLY("./F104/2009 .csv", percentage= 100)
F114_mcp.08<-mcp_analysis.POLY("./F114/2008 .csv", percentage= 100)
F114_mcp.09<-mcp_analysis.POLY("./F114/2009 .csv", percentage= 100)
F114_mcp.10<-mcp_analysis.POLY("./F114/2010 .csv", percentage= 100)
F114_mcp.11<-mcp_analysis.POLY("./F114/2011 .csv", percentage= 100)
F114_mcp.12<-mcp_analysis.POLY("./F114/2012 .csv", percentage= 100)
F137_mcp.09<-mcp_analysis.POLY("./F137/2009 .csv", percentage= 100)
F137_mcp.10<-mcp_analysis.POLY("./F137/2010 .csv", percentage= 100)
F137_mcp.11<-mcp_analysis.POLY("./F137/2011 .csv", percentage= 100)
F147_mcp.09<-mcp_analysis.POLY("./F147/2009 .csv", percentage= 100)
F147_mcp.10<-mcp_analysis.POLY("./F147/2010 .csv", percentage= 100)
F147_mcp.11<-mcp_analysis.POLY("./F147/2011 .csv", percentage= 100)
F147_mcp.12<-mcp_analysis.POLY("./F147/2012 .csv", percentage= 100)
F36_mcp.08<-mcp_analysis.POLY("./F36/2008 .csv", percentage= 100)
F36_mcp.09<-mcp_analysis.POLY("./F36/2009 .csv", percentage= 100)
F36_mcp.10<-mcp_analysis.POLY("./F36/2010 .csv", percentage= 100)
F36_mcp.11<-mcp_analysis.POLY("./F36/2011 .csv", percentage= 100)
F36_mcp.12<-mcp_analysis.POLY("./F36/2012 .csv", percentage= 100)
F66_mcp.08<-mcp_analysis.POLY("./F66/2008 .csv", percentage= 100)
F66_mcp.09<-mcp_analysis.POLY("./F66/2009 .csv", percentage= 100)
F66_mcp.10<-mcp_analysis.POLY("./F66/2010 .csv", percentage= 100)
M119_mcp.08<-mcp_analysis.POLY("./M119/2008 .csv", percentage= 100)
M119_mcp.09<-mcp_analysis.POLY("./M119/2009 .csv", percentage= 100)
M119_mcp.10<-mcp_analysis.POLY("./M119/2010 .csv", percentage= 100)
M112_mcp.07<-mcp_analysis.POLY("./M112/2007 .csv", percentage= 100)
M112_mcp.09<-mcp_analysis.POLY("./M112/2009 .csv", percentage= 100)
M112_mcp.10<-mcp_analysis.POLY("./M112/2010 .csv", percentage= 100)
M69_mcp.09<-mcp_analysis.POLY("./M69/2009 .csv", percentage= 100)
M69_mcp.10<-mcp_analysis.POLY("./M69/2010 .csv", percentage= 100)

# vignette("ggplot2-specs")

## Fortify mcp polygons for ggplot2 *YEAR*:
F104_mcp.08T <- fortify(F104_mcp.08, region = "id")
F104_mcp.09T <- fortify(F104_mcp.09, region = "id")
F114_mcp.08T <- fortify(F114_mcp.08, region = "id")
F114_mcp.09T <- fortify(F114_mcp.09, region = "id")
F114_mcp.10T <- fortify(F114_mcp.10, region = "id")
F114_mcp.11T <- fortify(F114_mcp.11, region = "id")
F114_mcp.12T <- fortify(F114_mcp.12, region = "id")
F137_mcp.09T <- fortify(F137_mcp.09, region = "id")
F137_mcp.10T <- fortify(F137_mcp.10, region = "id")
F137_mcp.11T <- fortify(F137_mcp.11, region = "id")
F147_mcp.09T <- fortify(F147_mcp.09, region = "id")
F147_mcp.10T <- fortify(F147_mcp.10, region = "id")
F147_mcp.11T <- fortify(F147_mcp.11, region = "id")
F147_mcp.12T <- fortify(F147_mcp.12, region = "id")
F36_mcp.08T <- fortify(F36_mcp.08, region = "id")
F36_mcp.09T <- fortify(F36_mcp.09, region = "id")
F36_mcp.10T <- fortify(F36_mcp.10, region = "id")
F36_mcp.11T <- fortify(F36_mcp.11, region = "id")
F36_mcp.12T <- fortify(F36_mcp.12, region = "id")
F66_mcp.08T <- fortify(F66_mcp.08, region = "id")
F66_mcp.09T <- fortify(F66_mcp.09, region = "id")
F66_mcp.10T <- fortify(F66_mcp.10, region = "id")
M119_mcp.08T <- fortify(M119_mcp.08, region = "id")
M119_mcp.09T <- fortify(M119_mcp.09, region = "id")
M119_mcp.10T <- fortify(M119_mcp.10, region = "id")
M112_mcp.07T <- fortify(M112_mcp.07, region = "id")
M112_mcp.09T <- fortify(M112_mcp.09, region = "id")
M112_mcp.10T <- fortify(M112_mcp.10, region = "id")
M69_mcp.09T <- fortify(M69_mcp.09, region = "id")
M69_mcp.10T <- fortify(M69_mcp.10, region = "id")
M215_mcp.11T <- fortify(M215_mcp.11, region = "id")
M215_mcp.12T <- fortify(M215_mcp.12, region = "id")


## Plot individual HR shift:
# mcp.shift.TEST3 <- ggplot() +
#   geom_polygon(data=F104_mcp.08T, aes(x=F104_mcp.08T$long, y=F104_mcp.08T$lat),
#                alpha=0.2,colour="black") +
#   geom_polygon(data=F104_mcp.09T, aes(x=F104_mcp.09T$long, y=F104_mcp.09T$lat),
#                alpha=0.2,colour="black") +
#   theme_bw() +labs(x="Easting (m)", y="Northing (m)",title="F104 Home Range Shift") +
#   theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5))
# mcp.shift.TEST3

mcp.shift.TEST4 <- ggplot() +
  geom_polygon(data=F104_mcp.08T, aes(x=F104_mcp.08T$long, y=F104_mcp.08T$lat),
               alpha=0.1,colour="black",linetype=2) +
  geom_polygon(data=F104_mcp.09T, aes(x=F104_mcp.09T$long, y=F104_mcp.09T$lat),
               alpha=0.1,colour="black",linetype=2) +
  geom_polygon(data=F114_mcp.08T, aes(x=F114_mcp.08T$long, y=F114_mcp.08T$lat),
               alpha=0.1,colour="black",linetype=3) +
  geom_polygon(data=F114_mcp.09T, aes(x=F114_mcp.09T$long, y=F114_mcp.09T$lat),
               alpha=0.1,colour="black",linetype=3) +
  geom_polygon(data=F114_mcp.10T, aes(x=F114_mcp.10T$long, y=F114_mcp.10T$lat),
               alpha=0.1,colour="black",linetype=3) +
  geom_polygon(data=F114_mcp.11T, aes(x=F114_mcp.11T$long, y=F114_mcp.11T$lat),
               alpha=0.1,colour="black",linetype=3) +
  geom_polygon(data=F114_mcp.12T, aes(x=F114_mcp.12T$long, y=F114_mcp.12T$lat),
               alpha=0.1,colour="black",linetype=3) +
  geom_polygon(data=F137_mcp.09T, aes(x=F137_mcp.09T$long, y=F137_mcp.09T$lat),
               alpha=0.1,colour="black",linetype=4) +
  geom_polygon(data=F137_mcp.10T, aes(x=F137_mcp.10T$long, y=F137_mcp.10T$lat),
               alpha=0.1,colour="black",linetype=4) +
  geom_polygon(data=F137_mcp.11T, aes(x=F137_mcp.11T$long, y=F137_mcp.11T$lat),
               alpha=0.1,colour="black",linetype=4) +
  geom_polygon(data=F147_mcp.09T, aes(x=F147_mcp.09T$long, y=F147_mcp.09T$lat),
               alpha=0.1,colour="red",linetype=1) +
  geom_polygon(data=F147_mcp.10T, aes(x=F147_mcp.10T$long, y=F147_mcp.10T$lat),
               alpha=0.1,colour="red",linetype=1) +
  geom_polygon(data=F147_mcp.11T, aes(x=F147_mcp.11T$long, y=F147_mcp.11T$lat),
               alpha=0.1,colour="red",linetype=1) +
  geom_polygon(data=F147_mcp.12T, aes(x=F147_mcp.12T$long, y=F147_mcp.12T$lat),
               alpha=0.1,colour="red",linetype=1) +
  geom_polygon(data=F36_mcp.08T, aes(x=F36_mcp.08T$long, y=F36_mcp.08T$lat),
               alpha=0.1,colour="black",linetype=6) +
  geom_polygon(data=F36_mcp.09T, aes(x=F36_mcp.09T$long, y=F36_mcp.09T$lat),
               alpha=0.1,colour="black",linetype=6) +
  geom_polygon(data=F36_mcp.10T, aes(x=F36_mcp.10T$long, y=F36_mcp.10T$lat),
               alpha=0.1,colour="black",linetype=6) +
  geom_polygon(data=F36_mcp.11T, aes(x=F36_mcp.11T$long, y=F36_mcp.11T$lat),
               alpha=0.1,colour="black",linetype=6) +
  geom_polygon(data=F36_mcp.12T, aes(x=F36_mcp.12T$long, y=F36_mcp.12T$lat),
               alpha=0.1,colour="black",linetype=6) +
  geom_polygon(data=F66_mcp.08T, aes(x=F66_mcp.08T$long, y=F66_mcp.08T$lat),
               alpha=0.1,colour="black",linetype=1) +
  geom_polygon(data=F66_mcp.09T, aes(x=F66_mcp.09T$long, y=F66_mcp.09T$lat),
               alpha=0.1,colour="black",linetype=1) +
  geom_polygon(data=F66_mcp.10T, aes(x=F66_mcp.10T$long, y=F66_mcp.10T$lat),
               alpha=0.1,colour="black",linetype=1) +
  geom_polygon(data=M119_mcp.08T, aes(x=M119_mcp.08T$long, y=M119_mcp.08T$lat),
               alpha=0.1,colour="blue",linetype=2) +
  geom_polygon(data=M119_mcp.09T, aes(x=M119_mcp.09T$long, y=M119_mcp.09T$lat),
               alpha=0.1,colour="blue",linetype=2) +
  geom_polygon(data=M119_mcp.10T, aes(x=M119_mcp.10T$long, y=M119_mcp.10T$lat),
               alpha=0.1,colour="blue",linetype=2) +
  geom_polygon(data=M112_mcp.07T, aes(x=M112_mcp.07T$long, y=M112_mcp.07T$lat),
               alpha=0.1,colour="blue",linetype=3) +
  geom_polygon(data=M112_mcp.09T, aes(x=M112_mcp.09T$long, y=M112_mcp.09T$lat),
               alpha=0.1,colour="blue",linetype=3) +
  geom_polygon(data=M112_mcp.10T, aes(x=M112_mcp.10T$long, y=M112_mcp.10T$lat),
               alpha=0.1,colour="blue",linetype=3) +
  # geom_polygon(data=M69_mcp.09T, aes(x=M69_mcp.09T$long, y=M69_mcp.09T$lat),
  #              alpha=0.1,colour="black") +
  # geom_polygon(data=M69_mcp.10T, aes(x=M69_mcp.10T$long, y=M69_mcp.10T$lat),
  #              alpha=0.1,colour="black") +
  # geom_polygon(data=M215_mcp.11T, aes(x=M215_mcp.11T$long, y=M215_mcp.11T$lat),
  #              alpha=0.1,colour="black") +
  # geom_polygon(data=M215_mcp.12T, aes(x=M215_mcp.12T$long, y=M215_mcp.12T$lat),
  #              alpha=0.1,colour="black") +
  theme_bw() +labs(x="Easting (m)", y="Northing (m)",title="Yearly Home Range Shifts") +
  theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5))

mcp.shift.TEST4
?geom_polygon
?ggplot



###########################################################################################
##                                       SEASON                                          ##


## Create MCP polygons by SEASON:
F66_mcp.EM<-mcp_analysis.POLY("./Emergence .csv", percentage= 100)
F66_mcp.DRY<-mcp_analysis.POLY("./Dry .csv", percentage= 100)
F66_mcp.MON<-mcp_analysis.POLY("./Monsoon .csv", percentage= 100)
F66_mcp.PM<-mcp_analysis.POLY("./Post_Monsoon .csv", percentage= 100)


## Fortify mcp polygons for ggplot2 *SEASON*:
M215_mcp.EMT <- fortify(M215_mcp.EM, region = "id")
M215_mcp.DRYT <- fortify(M215_mcp.DRY, region = "id")
M215_mcp.MONT <- fortify(M215_mcp.MON, region = "id")
M215_mcp.PMT <- fortify(M215_mcp.PM, region = "id")

M112_mcp.DRYT <- fortify(M112_mcp.DRY, region = "id")
M112_mcp.MONT <- fortify(M112_mcp.MON, region = "id")
M112_mcp.PMT <- fortify(M112_mcp.PM, region = "id")

M119_mcp.DRYT <- fortify(M119_mcp.DRY, region = "id")
M119_mcp.MONT <- fortify(M119_mcp.MON, region = "id")
M119_mcp.PMT <- fortify(M119_mcp.PM, region = "id")

F114_mcp.EMT <- fortify(F114_mcp.EM, region = "id")
F114_mcp.DRYT <- fortify(F114_mcp.DRY, region = "id")
F114_mcp.MONT <- fortify(F114_mcp.MON, region = "id")
F114_mcp.PMT <- fortify(F114_mcp.PM, region = "id")

F137_mcp.EMT <- fortify(F137_mcp.EM, region = "id")
F137_mcp.DRYT <- fortify(F137_mcp.DRY, region = "id")
F137_mcp.MONT <- fortify(F137_mcp.MON, region = "id")
F137_mcp.PMT <- fortify(F137_mcp.PM, region = "id")

F147_mcp.EMT <- fortify(F147_mcp.EM, region = "id")
F147_mcp.DRYT <- fortify(F147_mcp.DRY, region = "id")
F147_mcp.MONT <- fortify(F147_mcp.MON, region = "id")
F147_mcp.PMT <- fortify(F147_mcp.PM, region = "id")

F252_mcp.EMT <- fortify(F252_mcp.EM, region = "id")
F252_mcp.DRYT <- fortify(F252_mcp.DRY, region = "id")
F252_mcp.MONT <- fortify(F252_mcp.MON, region = "id")
F252_mcp.PMT <- fortify(F252_mcp.PM, region = "id")

F36_mcp.EMT <- fortify(F36_mcp.EM, region = "id")
F36_mcp.DRYT <- fortify(F36_mcp.DRY, region = "id")
F36_mcp.MONT <- fortify(F36_mcp.MON, region = "id")
F36_mcp.PMT <- fortify(F36_mcp.PM, region = "id")

F66_mcp.EMT <- fortify(F66_mcp.EM, region = "id")
F66_mcp.DRYT <- fortify(F66_mcp.DRY, region = "id")
F66_mcp.MONT <- fortify(F66_mcp.MON, region = "id")
F66_mcp.PMT <- fortify(F66_mcp.PM, region = "id")

mcp.shift.TEST5 <- ggplot() +
  geom_polygon(data=F114_mcp.EMT, aes(x=F114_mcp.EMT$long, y=F114_mcp.EMT$lat),
               alpha=0.1,colour="blue",linetype=2) +
  geom_polygon(data=F114_mcp.DRYT, aes(x=F114_mcp.DRYT$long, y=F114_mcp.DRYT$lat),
               alpha=0.1,colour="red",linetype=3) +
  geom_polygon(data=F114_mcp.MONT, aes(x=F114_mcp.MONT$long, y=F114_mcp.MONT$lat),
               alpha=0.1,colour="green",linetype=4) +
  geom_polygon(data=F114_mcp.PMT, aes(x=F114_mcp.PMT$long, y=F114_mcp.PMT$lat),
               alpha=0.1,colour="black",linetype=5) +
  geom_polygon(data=F137_mcp.EMT, aes(x=F137_mcp.EMT$long, y=F137_mcp.EMT$lat),
               alpha=0.1,colour="blue",linetype=2) +
  geom_polygon(data=F137_mcp.DRYT, aes(x=F137_mcp.DRYT$long, y=F137_mcp.DRYT$lat),
               alpha=0.1,colour="red",linetype=3) +
  geom_polygon(data=F137_mcp.MONT, aes(x=F137_mcp.MONT$long, y=F137_mcp.MONT$lat),
               alpha=0.1,colour="green",linetype=4) +
  geom_polygon(data=F137_mcp.PMT, aes(x=F137_mcp.PMT$long, y=F137_mcp.PMT$lat),
               alpha=0.1,colour="black",linetype=5) +
  geom_polygon(data=F147_mcp.EMT, aes(x=F147_mcp.EMT$long, y=F147_mcp.EMT$lat),
               alpha=0.1,colour="blue",linetype=2) +
  geom_polygon(data=F147_mcp.DRYT, aes(x=F147_mcp.DRYT$long, y=F147_mcp.DRYT$lat),
               alpha=0.1,colour="red",linetype=3) +
  geom_polygon(data=F147_mcp.MONT, aes(x=F147_mcp.MONT$long, y=F147_mcp.MONT$lat),
               alpha=0.1,colour="green",linetype=4) +
  geom_polygon(data=F147_mcp.PMT, aes(x=F147_mcp.PMT$long, y=F147_mcp.PMT$lat),
               alpha=0.1,colour="black",linetype=5) +
  # geom_polygon(data=F252_mcp.EMT, aes(x=F252_mcp.EMT$long, y=F252_mcp.EMT$lat),
  #              alpha=0.1,colour="black",linetype=2) +
  # geom_polygon(data=F252_mcp.DRYT, aes(x=F252_mcp.DRYT$long, y=F252_mcp.DRYT$lat),
  #              alpha=0.1,colour="black",linetype=3) +
  # geom_polygon(data=F252_mcp.MONT, aes(x=F252_mcp.MONT$long, y=F252_mcp.MONT$lat),
  #              alpha=0.1,colour="black",linetype=4) +
  # geom_polygon(data=F252_mcp.PMT, aes(x=F252_mcp.PMT$long, y=F252_mcp.PMT$lat),
  #              alpha=0.1,colour="black",linetype=5) +
  geom_polygon(data=F36_mcp.EMT, aes(x=F36_mcp.EMT$long, y=F36_mcp.EMT$lat),
               alpha=0.1,colour="blue",linetype=2) +
  geom_polygon(data=F36_mcp.DRYT, aes(x=F36_mcp.DRYT$long, y=F36_mcp.DRYT$lat),
               alpha=0.1,colour="red",linetype=3) +
  geom_polygon(data=F36_mcp.MONT, aes(x=F36_mcp.MONT$long, y=F36_mcp.MONT$lat),
               alpha=0.1,colour="green",linetype=4) +
  geom_polygon(data=F36_mcp.PMT, aes(x=F36_mcp.PMT$long, y=F36_mcp.PMT$lat),
               alpha=0.1,colour="black",linetype=5) +
  geom_polygon(data=F66_mcp.EMT, aes(x=F66_mcp.EMT$long, y=F66_mcp.EMT$lat),
               alpha=0.1,colour="blue",linetype=2) +
  geom_polygon(data=F66_mcp.DRYT, aes(x=F66_mcp.DRYT$long, y=F66_mcp.DRYT$lat),
               alpha=0.1,colour="red",linetype=3) +
  geom_polygon(data=F66_mcp.MONT, aes(x=F66_mcp.MONT$long, y=F66_mcp.MONT$lat),
               alpha=0.1,colour="green",linetype=4) +
  geom_polygon(data=F66_mcp.PMT, aes(x=F66_mcp.PMT$long, y=F66_mcp.PMT$lat),
               alpha=0.1,colour="black",linetype=5) +
  theme_bw() +labs(x="Easting (m)", y="Northing (m)",title="Seasonal Home Range Shifts") +
  theme(legend.position="none", plot.title = element_text(face = "bold", hjust = 0.5))

mcp.shift.TEST5 





#################################
##
##     SEASONAL REFUGE USE
##
#################################


#############################################
## REFUGE SPATIAL POINTS OF ALL GILA MONSTERS:
#############################################

# Plot spdf of refugia:
Refugia.df <- read.csv('./Refuge_Use/Refuge_Input_Red.csv')
View(Refugia.df)

# Using the function "SpatialPoints" we create an object of class SpatialPoints.
# We have to specify the coordinates, whereas the bbox is automatically generated.
Refugia.sp <- SpatialPoints(Refugia.df[,c("EASTING","NORTHING")], proj4string=CRS.SC)

## CREATE SPDF:
Refugia.spdf <- SpatialPointsDataFrame(coords=Refugia.df[,6:7],data=Refugia.df[,1:5],
                                      proj4string=CRS.SC)

Refigia_F <- fortify(Refugia.spdf, region = "SEASON")

## Plot Refuge Points:
Refugia.All <- ggplot() +
  geom_point(data=Refigia_F, aes(x=Refigia_F$long, y=Refigia_F$lat),
             shape=17)+
  theme_bw() +labs(x="Easting (m)", y="Northing (m)",title="Refuge Points of All Gila Monsters") +
  theme(legend.position=c(0.75, 0.8), plot.title = element_text(face = "bold", hjust = 0.5))

Refugia.All

## spplot() refuge sites:
spplot(Refugia.spdf, zcol="SEASON")

# custom color scales
spplot(Refugia.spdf, zcol="SEASON", col.regions=rainbow(4))
spplot(Refugia.spdf, zcol="SEASON", cex=.5, col.regions=rainbow(4))





#####################################################
##          Plot spatial points on raster image
##
#####################################################

# mcp_analysis.POLY <- function(filename, percentage){
  #   data <- read.csv(file = filename,stringsAsFactors = FALSE)
  #   data.sp <- data[, c("LIZARDNUMBER", "EASTING", "NORTHING")]
  #   coordinates(data.sp) <- c("EASTING", "NORTHING")
  #   proj4string(data.sp) <- CRS.SC
  #   mcp_out <- mcp(data.sp, percentage, unout="ha")
  # }

# F104_MCP<-mcp_analysis.POLY('./F104/F104 .csv', percentage= 100)

# F104.a <- read_csv("./F104/F104 .csv")
# F104_mcp.F <- fortify(F104_MCP, region = "id")
All.Gilas <- read_csv("./GM_Final_Data.csv")
View(All.Gilas)

# utm_points <- cbind(F104.a$EASTING, F104.a$NORTHING)
# F104.latlong <- cbind(F104_mcp.F$long, F104_mcp.F$lat)
# colnames(F104.latlong) <- c("x","y")
utm_points <- cbind(All.Gilas$EASTING, All.Gilas$NORTHING)

utm_locations <- SpatialPoints(utm_points, proj4string=CRS.SC)
# F104MCP_locations <- SpatialPolygons(F104.latlong, proj4string=CRS.SC)


proj_lat.lon <- as.data.frame(spTransform(utm_locations, CRS("+proj=longlat +datum=WGS84")))
colnames(proj_lat.lon) <- c("x","y")
# F104.MCP_lat.lon <- as.data.frame(spTransform(F104.latlong, CRS("+proj=longlat +datum=WGS84")))
# colnames(F104.MCP_lat.lon) <- c("x","y")

# raster_myMap<- ggmap(myMap, projection = CRS.SC)

## FORTIGY SPATIAL SPATIAL POINTS FOR PLOTTING:
proj_lat.lon <- fortify(proj_lat.lon, region = "Type")


# proj_lat.lon <- as.data.frame(spTransform(utm_locations, CRS("+proj=longlat +datum=WGS84")))
# colnames(proj_lat.lon) <- c("x","y")

## FORTIFY SPATIAL SPATIAL POLYGONS FOR PLOTTING:
# F104_mcp.F <- fortify(F104.MCP_lat.lon, region = "id")



myMap <- get_stamenmap(bbox = c(left = -111.009,
                                bottom = 32.459,
                                right = -110.969,
                                top = 32.474),
                       maptype = "terrain", 
                       crop = FALSE,
                       zoom = 15)


# plot map
# ggmap(myMap)+geom_point(data=proj_lat.lon, aes(x=x, y=y), size=1)
  # geom_polygon(data=F104_mcp.F, aes(x=F104_mcp.F$long, y=F104_mcp.F$lat),
  #              alpha=0.1,colour="blue",linetype=2)

ggmap(myMap)+geom_point(data=proj_lat.lon, aes(x=x, y=y), size=0.3)
# geom_polygon(data=F104_mcp.F, aes(x=F104_mcp.F$long, y=F104_mcp.F$lat),
#              alpha=0.1,colour="blue",linetype=2)

#####################################################
##          Plot MCPs on raster image
##
#####################################################

### attempt 1


# ggmap(APSU_SM) + 
#   geom_point(data=campus_points, aes(x = X, y = Y, color=Name), size = 4, alpha = 0.8) + 
#   geom_polygon(data = fortify(outline_poly), aes(long, lat, group = group), 
#                colour = "black", fill = NA, alpha = 0.5)

# rm(rasterTEST,utm_locations.TEST,proj_lat.lon.TEST,utm_points.TEST)

###

## Get/view the stamen map (bbox should be adjusted appropriately):
myMap <- get_stamenmap(bbox = c(left = -111.009,
                                bottom = 32.459,
                                right = -110.969,
                                top = 32.474),
                       maptype = "terrain", 
                       crop = FALSE,
                       zoom = 15)
 
myMap_imagery <- get_stamenmap(bbox = c(left = -111.009,
                                bottom = 32.459,
                                right = -110.969,
                                top = 32.474),
                       maptype = "toner-2011", 
                       crop = FALSE,
                       zoom = 15)
?get_stamenmap
ggmap(myMap_imagery)

Season.Map <- get_stamenmap(bbox = c(left = -111.005,
                                bottom = 32.46,
                                right = -110.98,
                                top = 32.475),
                       maptype = "toner-background", 
                       crop = FALSE,
                       zoom = 16)

                          
ggmap(Season.Map)
Season.Map2 <- get_stamenmap(bbox = c(left = -111.005,
                                      bottom = 32.46,
                                      right = -110.98,
                                      top = 32.475),
                             maptype = "terrain", 
                             crop = FALSE,
                             zoom = 16)

get_stamen_tile_download_fail_log()
ggmap(Season.Map2)
?get_stamenmap
rm(Season.Map2)
#toner-background

## The MCP I created had the easting/northing values but didn’t have the projection set 
## (see: mcp.out@proj4string, where mcp.out is the name of your MCP object for any given 
## animal). So first I set the polygon projection with proj4string() and then reprojected 
## the polygon to lat/lon with spTransform():
  
# F104_MCP@proj4string
# proj4string(mcp.out) <- CRS("+proj=utm +zone=12 +datum=WGS84")
F104_latlon <- spTransform(F104_MCP, CRS("+proj=longlat +datum=WGS84"))
F114_latlon <- spTransform(F114_MCP, CRS("+proj=longlat +datum=WGS84"))
F135_latlon <- spTransform(F135_MCP, CRS("+proj=longlat +datum=WGS84"))
F137_latlon <- spTransform(F137_MCP, CRS("+proj=longlat +datum=WGS84"))
F146_latlon <- spTransform(F146_MCP, CRS("+proj=longlat +datum=WGS84"))
F147_latlon <- spTransform(F147_MCP, CRS("+proj=longlat +datum=WGS84"))
F200_latlon <- spTransform(F200_MCP, CRS("+proj=longlat +datum=WGS84"))
F214_latlon <- spTransform(F214_MCP, CRS("+proj=longlat +datum=WGS84"))
F252_latlon <- spTransform(F252_MCP, CRS("+proj=longlat +datum=WGS84"))
F36_latlon <- spTransform(F36_MCP, CRS("+proj=longlat +datum=WGS84"))
F66_latlon <- spTransform(F66_MCP, CRS("+proj=longlat +datum=WGS84"))
M112_latlon <- spTransform(M112_MCP, CRS("+proj=longlat +datum=WGS84"))
M119_latlon <- spTransform(M119_MCP, CRS("+proj=longlat +datum=WGS84"))
M14_latlon <- spTransform(M14_MCP, CRS("+proj=longlat +datum=WGS84"))
M215_latlon <- spTransform(M215_MCP, CRS("+proj=longlat +datum=WGS84"))
M255_latlon <- spTransform(M255_MCP, CRS("+proj=longlat +datum=WGS84"))
M67_latlon <- spTransform(M67_MCP, CRS("+proj=longlat +datum=WGS84"))
M69_latlon <- spTransform(M69_MCP, CRS("+proj=longlat +datum=WGS84"))


## If your polygon has a projection then you can skip that first step. This gives you a 
## useable stamen map and a MCP polygon in lat/lon. Then all you need to do is use ggmap() 
## to map them:
  
SC_stamen_map <- ggmap(myMap) +
  # geom_point(data = proj_lat.lon, aes(x=x, y=y), size = 0.3, alpha = 0.8, color = "black") +
  geom_polygon(data = fortify(F104_latlon), aes(long, lat, group=group), colour = "red", 
             fill = NA) +
  geom_polygon(data = fortify(F114_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F135_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F137_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F146_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F147_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F200_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F214_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F252_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F36_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F66_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(M112_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M119_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M14_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M215_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M255_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M67_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M69_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  xlab("Longitude")+
  ylab("Latitude")
    
SC_stamen_map

SC_stamen_mapTONER <- ggmap(myMap_imagery) +
  # geom_point(data = proj_lat.lon, aes(x=x, y=y), size = 0.3, alpha = 0.8, color = "black") +
  geom_polygon(data = fortify(F104_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F114_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F135_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F137_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F146_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F147_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F200_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F214_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F252_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F36_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F66_latlon), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(M112_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M119_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M14_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M215_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M255_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M67_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M69_latlon), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  xlab("Longitude")+
  ylab("Latitude")

SC_stamen_mapTONER

F104_latlonK <- spTransform(F104_KDE, CRS("+proj=longlat +datum=WGS84"))
F114_latlonK <- spTransform(F114_KDE, CRS("+proj=longlat +datum=WGS84"))
F135_latlonK <- spTransform(F135_KDE, CRS("+proj=longlat +datum=WGS84"))
F137_latlonK <- spTransform(F137_KDE, CRS("+proj=longlat +datum=WGS84"))
F146_latlonK <- spTransform(F146_KDE, CRS("+proj=longlat +datum=WGS84"))
F147_latlonK <- spTransform(F147_KDE, CRS("+proj=longlat +datum=WGS84"))
F200_latlonK <- spTransform(F200_KDE, CRS("+proj=longlat +datum=WGS84"))
F214_latlonK <- spTransform(F214_KDE, CRS("+proj=longlat +datum=WGS84"))
F252_latlonK <- spTransform(F252_KDE, CRS("+proj=longlat +datum=WGS84"))
F36_latlonK <- spTransform(F36_KDE, CRS("+proj=longlat +datum=WGS84"))
F66_latlonK <- spTransform(F66_KDE, CRS("+proj=longlat +datum=WGS84"))
M112_latlonK <- spTransform(M112_KDE, CRS("+proj=longlat +datum=WGS84"))
M119_latlonK <- spTransform(M119_KDE, CRS("+proj=longlat +datum=WGS84"))
M14_latlonK <- spTransform(M14_KDE, CRS("+proj=longlat +datum=WGS84"))
M215_latlonK <- spTransform(M215_KDE, CRS("+proj=longlat +datum=WGS84"))
M255_latlonK <- spTransform(M255_KDE, CRS("+proj=longlat +datum=WGS84"))
M67_latlonK <- spTransform(M67_KDE, CRS("+proj=longlat +datum=WGS84"))
M69_latlonK <- spTransform(M69_KDE, CRS("+proj=longlat +datum=WGS84"))

SC_stamen_mapK <- ggmap(myMap) +
  # geom_point(data = proj_lat.lon, aes(x=x, y=y), size = 0.3, alpha = 0.8, color = "black") +
  geom_polygon(data = fortify(F104_latlonK), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F114_latlonK), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F135_latlonK), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F137_latlonK), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F146_latlonK), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F147_latlonK), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F200_latlonK), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F214_latlonK), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F252_latlonK), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F36_latlonK), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(F66_latlonK), aes(long, lat, group=group), colour = "red", 
               fill = NA) +
  geom_polygon(data = fortify(M112_latlonK), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M119_latlonK), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M14_latlonK), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M215_latlonK), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M255_latlonK), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M67_latlonK), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  geom_polygon(data = fortify(M69_latlonK), aes(long, lat, group=group), colour = "blue", 
               fill = NA) +
  xlab("Longitude")+
  ylab("Latitude")

SC_stamen_mapK

# Seasonal HR MCP Map:

F114_mcp.EMS<-mcp_analysis.POLY('./F114/Emergence .csv', percentage= 100)
F114_mcp.DRYS<-mcp_analysis.POLY('./F114/Dry .csv', percentage= 100)
F114_mcp.MONS<-mcp_analysis.POLY('./F114/Monsoon .csv', percentage= 100)
F114_mcp.PMS<-mcp_analysis.POLY('./F114/Post_Monsoon .csv', percentage= 100)
F137_mcp.EMS<-mcp_analysis.POLY('./F137/Emergence .csv', percentage= 100)
F137_mcp.DRYS<-mcp_analysis.POLY('./F137/Dry .csv', percentage= 100)
F137_mcp.MONS<-mcp_analysis.POLY('./F137/Monsoon .csv', percentage= 100)
F137_mcp.PMS<-mcp_analysis.POLY('./F137/Post_Monsoon .csv', percentage= 100)
F147_mcp.EMS<-mcp_analysis.POLY('./F147/Emergence .csv', percentage= 100)
F147_mcp.DRYS<-mcp_analysis.POLY('./F147/Dry .csv', percentage= 100)
F147_mcp.MONS<-mcp_analysis.POLY('./F147/Monsoon .csv', percentage= 100)
F147_mcp.PMS<-mcp_analysis.POLY('./F147/Post_Monsoon .csv', percentage= 100)
F36_mcp.EMS<-mcp_analysis.POLY('./F36/Emergence .csv', percentage= 100)
F36_mcp.DRYS<-mcp_analysis.POLY('./F36/Dry .csv', percentage= 100)
F36_mcp.MONS<-mcp_analysis.POLY('./F36/Monsoon .csv', percentage= 100)
F36_mcp.PMS<-mcp_analysis.POLY('./F36/Post_Monsoon .csv', percentage= 100)
F66_mcp.EMS<-mcp_analysis.POLY('./F66/Emergence .csv', percentage= 100)
F66_mcp.DRYS<-mcp_analysis.POLY('./F66/Dry .csv', percentage= 100)
F66_mcp.MONS<-mcp_analysis.POLY('./F66/Monsoon .csv', percentage= 100)
F66_mcp.PMS<-mcp_analysis.POLY('./F66/Post_Monsoon .csv', percentage= 100)

F114_latlonE <- spTransform(F114_mcp.EMS, CRS("+proj=longlat +datum=WGS84"))
F114_latlonD <- spTransform(F114_mcp.DRYS, CRS("+proj=longlat +datum=WGS84"))
F114_latlonM <- spTransform(F114_mcp.MONS, CRS("+proj=longlat +datum=WGS84"))
F114_latlonP <- spTransform(F114_mcp.PMS, CRS("+proj=longlat +datum=WGS84"))
F137_latlonE <- spTransform(F137_mcp.EMS, CRS("+proj=longlat +datum=WGS84"))
F137_latlonD <- spTransform(F137_mcp.DRYS, CRS("+proj=longlat +datum=WGS84"))
F137_latlonM <- spTransform(F137_mcp.MONS, CRS("+proj=longlat +datum=WGS84"))
F137_latlonP <- spTransform(F137_mcp.PMS, CRS("+proj=longlat +datum=WGS84"))
F147_latlonE <- spTransform(F147_mcp.EMS, CRS("+proj=longlat +datum=WGS84"))
F147_latlonD <- spTransform(F147_mcp.DRYS, CRS("+proj=longlat +datum=WGS84"))
F147_latlonM <- spTransform(F147_mcp.MONS, CRS("+proj=longlat +datum=WGS84"))
F147_latlonP <- spTransform(F147_mcp.PMS, CRS("+proj=longlat +datum=WGS84"))
F36_latlonE <- spTransform(F36_mcp.EMS, CRS("+proj=longlat +datum=WGS84"))
F36_latlonD <- spTransform(F36_mcp.DRYS, CRS("+proj=longlat +datum=WGS84"))
F36_latlonM <- spTransform(F36_mcp.MONS, CRS("+proj=longlat +datum=WGS84"))
F36_latlonP <- spTransform(F36_mcp.PMS, CRS("+proj=longlat +datum=WGS84"))
F66_latlonE <- spTransform(F66_mcp.EMS, CRS("+proj=longlat +datum=WGS84"))
F66_latlonD <- spTransform(F66_mcp.DRYS, CRS("+proj=longlat +datum=WGS84"))
F66_latlonM <- spTransform(F66_mcp.MONS, CRS("+proj=longlat +datum=WGS84"))
F66_latlonP <- spTransform(F66_mcp.PMS, CRS("+proj=longlat +datum=WGS84"))

SC_stamen_mapS <- ggmap(Season.Map) +
  # geom_point(data = proj_lat.lon, aes(x=x, y=y), size = 0.3, alpha = 0.8, color = "black") +
  geom_polygon(data = fortify(F114_latlonE), aes(long, lat, group=group), colour = "black", 
               linetype=2, fill = NA) +
  geom_polygon(data = fortify(F114_latlonD), aes(long, lat, group=group), colour = "red", 
               linetype=2, fill = NA) +
  geom_polygon(data = fortify(F114_latlonM), aes(long, lat, group=group), colour = "blue", 
               linetype=2, fill = NA) +
  geom_polygon(data = fortify(F114_latlonP), aes(long, lat, group=group), colour = "green", 
               linetype=2, fill = NA) +
  geom_polygon(data = fortify(F137_latlonE), aes(long, lat, group=group), colour = "black", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F137_latlonD), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F137_latlonM), aes(long, lat, group=group), colour = "blue", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F137_latlonP), aes(long, lat, group=group), colour = "green", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F147_latlonE), aes(long, lat, group=group), colour = "black", 
               linetype=3, fill = NA) +
  geom_polygon(data = fortify(F147_latlonD), aes(long, lat, group=group), colour = "red", 
               linetype=3, fill = NA) +
  geom_polygon(data = fortify(F147_latlonM), aes(long, lat, group=group), colour = "blue", 
               linetype=3, fill = NA) +
  geom_polygon(data = fortify(F147_latlonP), aes(long, lat, group=group), colour = "green", 
               linetype=3, fill = NA) +
  geom_polygon(data = fortify(F36_latlonE), aes(long, lat, group=group), colour = "black", 
               linetype=4, fill = NA) +
  geom_polygon(data = fortify(F36_latlonD), aes(long, lat, group=group), colour = "red", 
               linetype=4, fill = NA) +
  geom_polygon(data = fortify(F36_latlonM), aes(long, lat, group=group), colour = "blue", 
               linetype=4, fill = NA) +
  geom_polygon(data = fortify(F36_latlonP), aes(long, lat, group=group), colour = "green", 
               linetype=4, fill = NA) +
  geom_polygon(data = fortify(F66_latlonE), aes(long, lat, group=group), colour = "black", 
               linetype=5, fill = NA) +
  geom_polygon(data = fortify(F66_latlonD), aes(long, lat, group=group), colour = "red", 
               linetype=5, fill = NA) +
  geom_polygon(data = fortify(F66_latlonM), aes(long, lat, group=group), colour = "blue", 
               linetype=5, fill = NA) +
  geom_polygon(data = fortify(F66_latlonP), aes(long, lat, group=group), colour = "green",
               linetype=5, fill = NA) +
  xlab("Longitude")+
  ylab("Latitude")

SC_stamen_mapS

## Yearly shift stamen maps: 

F114_mcp.08<-mcp_analysis.POLY("./F114/2008 .csv", percentage= 100)
F114_mcp.09<-mcp_analysis.POLY("./F114/2009 .csv", percentage= 100)
F114_mcp.10<-mcp_analysis.POLY("./F114/2010 .csv", percentage= 100)
F114_mcp.11<-mcp_analysis.POLY("./F114/2011 .csv", percentage= 100)
F114_mcp.12<-mcp_analysis.POLY("./F114/2012 .csv", percentage= 100)
F137_mcp.09<-mcp_analysis.POLY("./F137/2009 .csv", percentage= 100)
F137_mcp.10<-mcp_analysis.POLY("./F137/2010 .csv", percentage= 100)
F137_mcp.11<-mcp_analysis.POLY("./F137/2011 .csv", percentage= 100)
F147_mcp.09<-mcp_analysis.POLY("./F147/2009 .csv", percentage= 100)
F147_mcp.10<-mcp_analysis.POLY("./F147/2010 .csv", percentage= 100)
F147_mcp.11<-mcp_analysis.POLY("./F147/2011 .csv", percentage= 100)
F147_mcp.12<-mcp_analysis.POLY("./F147/2012 .csv", percentage= 100)
F66_mcp.08<-mcp_analysis.POLY("./F66/2008 .csv", percentage= 100)
F66_mcp.09<-mcp_analysis.POLY("./F66/2009 .csv", percentage= 100)
F66_mcp.10<-mcp_analysis.POLY("./F66/2010 .csv", percentage= 100)
M119_mcp.08<-mcp_analysis.POLY("./M119/2008 .csv", percentage= 100)
M119_mcp.09<-mcp_analysis.POLY("./M119/2009 .csv", percentage= 100)
M119_mcp.10<-mcp_analysis.POLY("./M119/2010 .csv", percentage= 100)
M112_mcp.07<-mcp_analysis.POLY("./M112/2007 .csv", percentage= 100)
M112_mcp.09<-mcp_analysis.POLY("./M112/2009 .csv", percentage= 100)
M112_mcp.10<-mcp_analysis.POLY("./M112/2010 .csv", percentage= 100)


F114_latlon.08 <- spTransform(F114_mcp.08, CRS("+proj=longlat +datum=WGS84"))
F114_latlon.09 <- spTransform(F114_mcp.09, CRS("+proj=longlat +datum=WGS84"))
F114_latlon.10 <- spTransform(F114_mcp.10, CRS("+proj=longlat +datum=WGS84"))
F114_latlon.11 <- spTransform(F114_mcp.11, CRS("+proj=longlat +datum=WGS84"))
F114_latlon.12 <- spTransform(F114_mcp.12, CRS("+proj=longlat +datum=WGS84"))
F137_latlon.09 <- spTransform(F137_mcp.09, CRS("+proj=longlat +datum=WGS84"))
F137_latlon.10 <- spTransform(F137_mcp.10, CRS("+proj=longlat +datum=WGS84"))
F137_latlon.11 <- spTransform(F137_mcp.11, CRS("+proj=longlat +datum=WGS84"))
F147_latlon.09 <- spTransform(F147_mcp.09, CRS("+proj=longlat +datum=WGS84"))
F147_latlon.10 <- spTransform(F147_mcp.10, CRS("+proj=longlat +datum=WGS84"))
F147_latlon.11 <- spTransform(F147_mcp.11, CRS("+proj=longlat +datum=WGS84"))
F147_latlon.12 <- spTransform(F147_mcp.12, CRS("+proj=longlat +datum=WGS84"))
F66_latlon.08 <- spTransform(F66_mcp.08, CRS("+proj=longlat +datum=WGS84"))
F66_latlon.09 <- spTransform(F66_mcp.09, CRS("+proj=longlat +datum=WGS84"))
F66_latlon.10 <- spTransform(F66_mcp.10, CRS("+proj=longlat +datum=WGS84"))
M119_latlon.08 <- spTransform(M119_mcp.08, CRS("+proj=longlat +datum=WGS84"))
M119_latlon.09 <- spTransform(M119_mcp.09, CRS("+proj=longlat +datum=WGS84"))
M119_latlon.10 <- spTransform(M119_mcp.10, CRS("+proj=longlat +datum=WGS84"))
M112_latlon.07 <- spTransform(M112_mcp.07, CRS("+proj=longlat +datum=WGS84"))
M112_latlon.09 <- spTransform(M112_mcp.09, CRS("+proj=longlat +datum=WGS84"))
M112_latlon.10 <- spTransform(M112_mcp.10, CRS("+proj=longlat +datum=WGS84"))

MCP.Shift.Yearly <- ggmap(Season.Map) +
  # geom_point(data = proj_lat.lon, aes(x=x, y=y), size = 0.3, alpha = 0.8, color = "black") +
  geom_polygon(data = fortify(F114_latlon.08), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F114_latlon.09), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F114_latlon.10), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F114_latlon.11), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F114_latlon.12), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  # geom_polygon(data = fortify(F137_latlon.09), aes(long, lat, group=group), colour = "red", 
  #              linetype=1, fill = NA) +
  # geom_polygon(data = fortify(F137_latlon.10), aes(long, lat, group=group), colour = "red", 
  #              linetype=1, fill = NA) +
  # geom_polygon(data = fortify(F137_latlon.11), aes(long, lat, group=group), colour = "red", 
  #              linetype=1, fill = NA) +
  geom_polygon(data = fortify(F147_latlon.09), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F147_latlon.10), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F147_latlon.11), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F147_latlon.12), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F66_latlon.08), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F66_latlon.09), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(F66_latlon.10), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(M119_latlon.08), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(M119_latlon.09), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  geom_polygon(data = fortify(M119_latlon.10), aes(long, lat, group=group), colour = "red", 
               linetype=1, fill = NA) +
  # geom_polygon(data = fortify(M112_latlon.07), aes(long, lat, group=group), colour = "red", 
  #              linetype=6, fill = NA) +
  # geom_polygon(data = fortify(M112_latlon.09), aes(long, lat, group=group), colour = "red",
  #              linetype=6, fill = NA) +
  # geom_polygon(data = fortify(M112_latlon.10), aes(long, lat, group=group), colour = "red",
  #              linetype=6, fill = NA) +
  xlab("Longitude") +
  ylab("Latitude")

MCP.Shift.Yearly

#####################################
## Add North arrow and scale bar to stamen map :

## add north arrow
library(ggsn)
# north(data = NULL, location = "topright", scale = 0.1, symbol = 1,
#       x.min, x.max, y.min, y.max, anchor = NULL)

# SC_stamen_map + north(data = NULL, location = "topright", scale = 0.1, symbol = 1,
#       0.65, -0.9, 0.5, 0.9, anchor = NULL)

north2(SC_stamen_map, x = 0.2, y = 0.27, scale = 0.1, symbol = 16)

## Add scale bare
# scalebar(data = NULL, location = "bottomright", dist = NULL,
#          dist_unit = NULL, transform = NULL, dd2km = NULL, model = NULL,
#          height = 0.02, st.dist = 0.02, st.bottom = TRUE, st.size = 5,
#          st.color = "black", box.fill = c("black", "white"),
#          box.color = "black", border.size = 1, x.min = NULL, x.max = NULL,
#          y.min = NULL, y.max = NULL, anchor = NULL, facet.var = NULL,
#          facet.lev = NULL, ...)
# SC_stamen_map +
#   scalebar(x.min = -110.01, x.max = -111.00,
#            y.min = 32.57, y.max = 32.60,
#            dist = 4, dist_unit = "km",
#            st.bottom = FALSE, st.color = "black",
#            transform = FALSE, model = "WGS84") +
#   north2(SC_stamen_map, x = 0.2, y = 0.27, scale = 0.1, symbol = 16)

# For data in UTM meters: 
 
#  Scalebar and north arrow for MCP Maps:  

SC_stamen_map<-SC_stamen_map + ggsn::scalebar(x.min = -110.972, x.max = -110.966,
                     y.min = 32.477, y.max = 32.479, 
                     dist = 500, 
                     dist_unit="m", 
                     height=0.19,
                     st.bottom=FALSE, 
                     st.dist=0.3,
                     st.size=3,
                     transform = TRUE, 
                     model = 'WGS84') 
SC_stamen_map

SC_stamen_map+north2(SC_stamen_map, x = 0.836, y = 0.73, scale = 0.1, symbol = 16)

# scalebar and north arrow for KDE Maps:

SC_stamen_mapK <- SC_stamen_mapK + ggsn::scalebar(x.min = -110.972, x.max = -110.966,
                                              y.min = 32.477, y.max = 32.479, 
                                              dist = 500, 
                                              dist_unit="m", 
                                              height=0.19,
                                              st.bottom=FALSE, 
                                              st.dist=0.3,
                                              st.size=3,
                                              transform = TRUE, 
                                              model = 'WGS84') 
SC_stamen_mapK

SC_stamen_mapK + north2(SC_stamen_mapK, x = 0.836, y = 0.73, scale = 0.1, symbol = 16)

# scalebar and north arrow for season Maps:

SC_stamen_mapS <- SC_stamen_mapS + ggsn::scalebar(x.min = -111.005, x.max = -110.998,
                                                  y.min = 32.461, y.max = 32.463, 
                                                  dist = 250, 
                                                  dist_unit="m", 
                                                  height=0.19,
                                                  st.bottom=FALSE, 
                                                  st.dist=0.3,
                                                  st.size=3,
                                                  transform = TRUE, 
                                                  model = 'WGS84') 
SC_stamen_mapS

SC_stamen_mapS + north2(SC_stamen_mapS, x = 0.36, y = 0.30, scale = 0.1, symbol = 16)

MCP.Shift.Yearly <- MCP.Shift.Yearly + ggsn::scalebar(x.min = -111.005, x.max = -110.998,
                                                  y.min = 32.461, y.max = 32.463, 
                                                  dist = 250, 
                                                  dist_unit="m", 
                                                  height=0.19,
                                                  st.bottom=FALSE, 
                                                  st.dist=0.3,
                                                  st.size=3,
                                                  transform = TRUE, 
                                                  model = 'WGS84') 
MCP.Shift.Yearly

MCP.Shift.Yearly + north2(MCP.Shift.Yearly, x = 0.36, y = 0.30, scale = 0.1, symbol = 16)

#######################################################
##
##    DISTRIBUTION MAP WITH GEOGRAPHIC ATTRIBUTES    ##
##
#######################################################
#Obtain FIPS Codes by county 
fips <- county.fips

#Create county polygons
florida <- map(database = "county", regions = "florida", fill=T, plot=F)
IDs <- sub("^florida,","",florida$names)

#Add FIPS codes to the county polygons
fips.codes <- separate(data = fips, col = polyname, into = c("state", "county"), sep = ",")
fl_fips <- subset(fips.codes, state=="florida", select=fips)
names <- fips.codes$county
fl_IDs <- unique(fl_fips$fips)

#Create spatial polygons
fl_sp = map2SpatialPolygons(florida,fl_fips$fips,CRS("+proj=longlat"))
names(fl_sp@polygons) <- fl_IDs

###############################
###############################

world <- map_data("world")
states <- map_data("state")
counties <- map_data("county")

counties$polyname <- paste(counties$region, counties$subregion, sep = ",")
counties <- counties %>% left_join(fips, by = c("polyname" = "polyname"))
counties$fips <- as.character(counties$fips)

southwestern_states <- subset(states, region %in% 
                            c("arizona", "california", "utah", "nevada", 
                              "new mexico", "colorado"))

southwestern_counties <- subset(counties, region %in% 
                              c("arizona", "california", "utah", "nevada", 
                                "new mexico", "colorado"))

# florida_counties <- subset(southern_counties, region == "florida")

dist_map <- ggplot() + 
  geom_polygon(data = world, aes(x=long,y=lat, group=group), fill = "gray95", color = "white") +
  geom_polygon(data = states, aes(x=long,y=lat, group=group), fill = "gray", color = "white") +
  # geom_polygon(data = fl_poly, aes(x=long, y=lat, group=group, fill = fill))  
  geom_polygon(data = southwestern_states, aes(x=long,y=lat, group=group), fill = NA, color = "white") +
  geom_polygon(data = southwestern_counties, aes(x=long,y=lat, group=group), fill = NA, color = "black", size = 0.05) +
  coord_map("conic", lat0 = 30, xlim=c(-120,-105), ylim=c(30,42)) +
  scale_fill_identity() +
  theme_grey() + theme(legend.position="right") + theme(legend.title.align=0.5) +
  theme(panel.background = element_rect(fill = 'deepskyblue'),
        panel.grid.major = element_line(colour = NA)) +
  labs(x = "Longitude", y = "Latitude", fill = "Child Poverty", 
       title = "U.S. Gila Monster Distribution") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
## use to preview the map
dist_map



################################################################################
##
##                             FIGURE FORGE
##
################################################################################

grid.arrange(Raw.YearHR, yr.mean.adj,Raw.KDE.HR, LSM.KDE.HR, nrow = 2)

library(ggpubr)
ggarrange(Raw.YearHR, yr.mean.adj,Raw.KDE.HR, LSM.KDE.HR, labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)

########################## SUB VS NONSUB KDEs ##################
View(year)
RMmod.95kde<-lmer(Home_Range_95kde~Environment+Year+Sex+N+Environment*Sex+(1|Gila),data = year)

RM.marginal <- lsmeans(RMmod.95kde, 
                       ~ Environment)
# RM.marginal

## CATAGORIZE LSM GRAPH BY SEX BETWEEN ENVIRONMENT:
refRM_kde <- lsmeans(RMmod.95kde, specs = c("Environment","Sex"))

# refRM_sex
ref_dfRM_kde <- as.data.frame(summary(refRM_kde))
pd_RM <- position_dodge(0.1)

kde.mean.adj<-ggplot(ref_dfRM_kde, aes(x=Sex, y=lsmean, group=Environment))+
  geom_point(aes(shape = factor(Environment)), size = 3, position=position_dodge(.1), 
             show.legend = FALSE)+
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.1, position=position_dodge())+
  xlab("")+
  ylab("")
kde.mean.adj

pd_RM <- position_dodge(0.1)

Raw.kde<-ggplot(YR_Means.95KDEall, aes(x=Sex, y=Home_Range_95kde, group=Environment))+
  geom_point(aes(shape = factor(Environment)), size = 3,position=position_dodge(.1),
             show.legend = FALSE)+
  geom_errorbar(aes(ymin=Home_Range_95kde-se, ymax=Home_Range_95kde+se),
                width=.1, position=position_dodge())+
  xlab("")+
  ylab("95% KDE Area (ha)")
Raw.kde
# Raw.kde<-Raw.kde + theme(legend.title = element_blank(),
#                                legend.text = element_text(size = 12),
#                                legend.justification=c(0,1),
#                                legend.position=c(0.05, 0.95),
#                                legend.background = element_blank(),
#                                legend.key = element_blank(),
#                                legend.box.background = element_rect(colour = "black")) +
#   scale_shape_discrete(name  ="",
#                        breaks=c("nonsubsidized", "subsidized"),
#                        labels=c("Nonsubsidized", "Subsidized"))
# Raw.kde

library(gridExtra)
library(grid)

ggarrange(Raw.YearHR, yr.mean.adj,Raw.kde, kde.mean.adj, labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)

ggarrange(raw.seasonal, adj.seasonal, labels = c("A", "B"),
          nrow = 2)

######################################
## 100% MCP Multiple reg:

Graph1<-ggplot(year,aes(x=N100,y=Home_Range_100mcp))+
  geom_point(aes(shape = factor(Environment)), size = 4)+
  scale_shape_manual(values=c(16, 2), name="", breaks=c("nonsubsidized", "subsidized"),
                     labels=c("Nonsubsidized", "Subsidized"))+
  geom_smooth(aes(linetype=Environment),colour="black", method="lm", show.legend =FALSE) +
  scale_linetype_manual(values=c("solid", "dashed"))+
  xlab("Number of Relocations")+
  ylab("100% MCP Area (ha)")+
  scale_x_continuous(limits= c(0,125), breaks = c(0,25,50,75,100,125)) +
  guides(color = guide_legend(override.aes = list(linetype = 0)))+
  theme(plot.caption = element_text(hjust = 0,lineheight = 0.9)) +
  theme_bw() +
  theme(legend.position = c(.87,.85), legend.background = element_rect(colour = "black"),
        axis.text.x=element_blank(),
        axis.text.y  = element_text(vjust=0.5, size=14),
        axis.title.y  = element_text(size=18),
        axis.title.x  = element_blank(),
        legend.text = element_text(size = 12, face = "bold"),
        strip.text = element_text(size=12)) 

# Graph1<-Graph1+theme(axis.title=element_text(size = 18))

# legend at top-left, inside the plot
SCOH.hr.fig<-Graph1 + theme(legend.title = element_blank(),
                            legend.text = element_text(size = 14),
                            legend.justification=c(0,1),
                            legend.position=c(0.05, 0.95),
                            legend.background = element_blank(),
                            legend.key = element_blank(),
                            legend.box.background = element_rect(colour = "black"))
  # scale_shape_discrete(name  ="",
  #                      breaks=c("nonsubsidized", "subsidized"),
  #                      labels=c("Nonsubsidized", "Subsidized")) +
  # scale_linetype_discrete(name  ="",
  #                          breaks=c("nonsubsidized", "subsidized"),
  #                          labels=c("Nonsubsidized", "Subsidized"))


SCOH.hr.fig

######################################
## 95% KDE Multiple reg:

Graph2<-ggplot(year2,aes(x=N,y=Home_Range_95kde))+
  geom_point(aes(shape = factor(Environment)), size = 4)+
  scale_shape_manual(values=c(16, 2), name="", breaks=c("nonsubsidized", "subsidized"),
                     labels=c("Nonsubsidized", "Subsidized"))+
  geom_smooth(aes(linetype=Environment),colour="black", method="lm", show.legend =FALSE) +
  scale_linetype_manual(values=c("solid", "dashed"))+
  xlab("Number of Relocations")+
  ylab("95% KDE Area (ha)")+
  guides(color = guide_legend(override.aes = list(linetype = 0)))+
  theme(plot.caption = element_text(hjust = 0,lineheight = 0.9)) +
  theme_bw() +
  theme(legend.position = "none", legend.background = element_rect(colour = "black"),
        axis.text.x=element_text(vjust=0.5, size=14),
        axis.text.y  = element_text(vjust=0.5, size=14),
        axis.title.y  = element_text(size=18),
        axis.title.x  = element_text(size=18),
        # legend.text = element_text(size = 12, face = "bold"),
        strip.text = element_text(size=12))

# Graph2<-Graph2+theme(axis.title=element_text(size = 18))

# legend at top-left, inside the plot
# SCOH.hr.fig2<-Graph2 + theme(legend.title = element_blank(),
#                             legend.text = element_text(size = 14),
#                             legend.justification=c(0,1),
#                             legend.position=c(0.05, 0.95),
#                             legend.background = element_blank(),
#                             legend.key = element_blank(),
#                             legend.box.background = element_rect(colour = "black"))
# scale_shape_discrete(name  ="",
#                      breaks=c("nonsubsidized", "subsidized"),
#                      labels=c("Nonsubsidized", "Subsidized")) +
# scale_linetype_discrete(name  ="",
#                          breaks=c("nonsubsidized", "subsidized"),
#                          labels=c("Nonsubsidized", "Subsidized"))


SCOH.hr.fig2<-Graph2 + 
  theme(legend.position = "none") + scale_x_continuous(limits= c(0,125), 
                                                       breaks = c(0,25,50,75,100,125))
SCOH.hr.fig2

ggarrange(SCOH.hr.fig, SCOH.hr.fig2, labels = c("A", "B"),
          ncol = 1)


#########################
## RECHECKING ASSUMPTIONS
View(year)
year2<-read_csv("GM_Consolidated_ByYear_Input.csv")
View(nonsub2)

library(lme4)
library(readr)
year <- read_csv("GM_Consolidated_ByYear.csv")


nonsub2 <- subset(year, Environment == "nonsubsidized")

NSub.yearAffect <- lmer(Home_Range_100mcp~Year+Sex+Year*Sex+(1|Gila), data=nonsub2)
summary(NSub.yearAffect)
ggqqplot(year)

library(rstatix)
ggqqplot(year,"Home_Range_100mcp",facet.by="Environment")

#########################################################
RMmod.year95<-lmer(Home_Range_95mcp~Year+N95+Environment*Sex+
                     (1|Gila),data = year)
summary(RMmod.year95)
anova(RMmod.year95)

RM.mod.Season <- lmer(Home_Range_100mcp~Environment+Season+Sex+N+Environment*Season+
                        Season*Sex+(1|Gila), data=seasonal)

summary(RM.mod.Season)
anova(RM.mod.Season)

## calc effect size
library(effsize)
cohen.d(Environment, hedges.correction=TRUE,data=year)


