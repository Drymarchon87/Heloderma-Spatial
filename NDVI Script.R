
#######################################################################
##
##                   NDVI for Gila Monster Study Sites
##                     Stone Canyon/Owl Head Buttes
##
#######################################################################

## Load required packages
library(raster)
library(rgdal)
library(ggplot2)

###############################################
## Load each band image to object for stacking
###############################################

# Blue
B2 <- raster('./Landsat/LC08_L1TP_036037_20160705_20170221_01_T1_B2.tif') 
# Green
B3 <- raster('./Landsat/LC08_L1TP_036037_20160705_20170221_01_T1_B3.tif')
# Red
B4 <- raster('./Landsat/LC08_L1TP_036037_20160705_20170221_01_T1_B4.tif')
# NIR
B5 <- raster('./Landsat/LC08_L1TP_036037_20160705_20170221_01_T1_B5.tif')
# # SWIR 1
# B6 <- raster('./Landsat/LC08_L1TP_036037_20160705_20170221_01_T1_B6.tif')
# # SWIR 2
# B7 <- raster('./Landsat/LC08_L1TP_036037_20160705_20170221_01_T1_B7.tif')
# # Panchromatic
# B8 <- raster('./Landsat/LC08_L1TP_036037_20160705_20170221_01_T1_B8.tif')
# # Cirrus
# B9 <- raster('./Landsat/LC08_L1TP_036037_20160705_20170221_01_T1_B8.tif')
# # Thermal Infrared (TIRS) 1
# B10 <- raster('./Landsat/LC08_L1TP_036037_20160705_20170221_01_T1_B10.tif')
# # Thermal Infrared (TIRS) 2
# B11 <- raster('./Landsat/LC08_L1TP_036037_20160705_20170221_01_T1_B11.tif')

# Create a raster stack object with bands used
# S <- stack(B5, B4, B3, B2)
# rm(s,S)

## Can also create the RasterStack and filenames: **DOESN'T WORK**
# filenames <- paste0('.NDVI/LC08_L1TP_036038_20130627_20180131_01_T1_B', 1:11, ".tif") 
# filenames
# landsat <- stack(filenames)
# landsat

# Create list of NDVI file paths **CAN TRY THIS INSTEAD OF THE ABOVE
# **NOT WORKING**
# all_HARV_NDVI <- list.files("NEON-DS-Landsat-NDVI/HARV/2011/NDVI",
#                             full.names = TRUE,
#                             pattern = ".tif$")



###############################################
## Create True and False Color Composites
###############################################

## Make a “true (or natural) color” image, something that looks like a normal
## photograph (vegetation in green, water blue etc), need bands in the red, green and blue 
## regions. For this Landsat image, band 4 (red), 3 (green), and 2 (blue). 
## The plotRGB method used to combine them into a single composite. Can also supply
## additional arguments to plotRGB to improve the visualization (e.g. a linear stretch of the 
## values, using strecth = "lin").

# True Color composite:
LANDSAT.RGB <- stack(B4, B3, B2)
plotRGB(LANDSAT.RGB, axes = TRUE, stretch = "lin", main = "Landsat True Color Composite")

## True-color composite reveals much more about the landscape than the earlier gray images. 
## Another image visualization method in remote sensing is a “false color” image, 
## NIR, red, and green bands are combined. Makes it easy to see vegetation (in red).

# False color composite:
# Make one row of two columns to visualize side by side...
par(mfrow = c(1,2))
# plot RGB (true color comp) first
plotRGB(LANDSAT.RGB, axes=TRUE, stretch="lin", main="Landsat True Color Composite")
# Create stack for false color comp using NIR, Red, and Green bands and plot
LANDSAT.FCC <- stack(B5, B4, B3)
plotRGB(LANDSAT.FCC, axes=TRUE, stretch="lin", main="Landsat False Color Composite")





###############################################
## Subset and Rename Bands
###############################################

# select first 3 bands only
# landsatsub1 <- subset(landsat, 1:3) # same
# landsatsub2 <- landsat[[1:3]]

## Can and useful to set the names of each band used
names(LANDSAT.FCC)
names(LANDSAT.FCC) <- c('NIR', 'Red', 'green') 
names(LANDSAT.FCC)
names(LANDSAT.RGB)
names(LANDSAT.RGB) <- c('Red', 'Green', 'Blue') 
names(LANDSAT.RGB)





##############################################
## SPATIAL SUBSETTING OR CROPPING:
##############################################

## Can use spatial subsetting to limit analysis to a geographic subset of the image. 
## Spatial subsets can be created with the crop function, using an extent object, or another 
## spatial object from which an Extent can be extracted.

# Check the extent of the image, will be the same for RGB or False color (FCC) from above
extent(LANDSAT.FCC)

###########################
# STONE CANYON landsat crop:

# assign the extent of the site for subsetting
E <- extent(499450, 502000, 3591500, 3593300)

# Crop the landsat by the assigned extent (cropped for both RGB and FCC)
landsatcrop.A <- crop(LANDSAT.RGB, E)
landsatcrop.B <- crop(LANDSAT.FCC, E)

# extent(landsatcropA) # cancheck the new extent
# reset the plots to 1x1
par(mfrow = c(1,1))
# plot the new subsetted RGB or FCC
plotRGB(landsatcrop.B, axes=TRUE, stretch="lin", main="Landsat False Color Composite")


##########################
# OWL HEAD landsat crop: ***** NOT WORKING *****

# assign the extent of the site for subsetting
E2 <- extent(487927, 490477, 3606200, 3608000)

# Crop the landsat by the assigned extent (cropped for both RGB and FCC)
landsatcrop.A2 <- crop(LANDSAT.RGB, E2)
landsatcrop.B2 <- crop(LANDSAT.FCC, E2)

# plot the new subsetted RGB or FCC
plotRGB(landsatcrop.B2, axes=TRUE, stretch="lin", main="Landsat False Color Composite")





##############################################
##      RELATIONSHIP BETWEEN BANDS
##############################################

## Scatterplot matrix can be helpful in exploring relationships between raster layers. 
## Can be done with the pairs() function of the raster package

## Plot of reflection in the ultra-blue wavelength against reflection in the blue wavelength.
# STONE CANYON
pairs(landsatcrop.A[[2:3]], main = "Blue versus Green")
pairs(landsatcrop.B[[1:2]], main = "Red versus NIR")
# OWL HEAD
pairs(landsatcrop.A2[[2:3]], main = "Blue versus Green")
pairs(landsatcrop.B2[[1:2]], main = "Red versus NIR")





##############################################
##                  NDVI
##############################################

# Create stacks ***** NOT WORKING *****
# raslist <- paste0('.NDVI/LC08_L1TP_036038_20130627_20180131_01_T1_B', 1:11, ".tif") 
# landsat <- stack(raslist)

# Manually create raster stacks if havnt aleady done
# Blue
b2 <- raster('./NDVI/LC08_L1TP_036038_20130627_20180131_01_T1_B2.tif') 
# Green
b3 <- raster('./NDVI/LC08_L1TP_036038_20130627_20180131_01_T1_B3.tif')
# Red
b4 <- raster('./NDVI/LC08_L1TP_036038_20130627_20180131_01_T1_B4.tif')
# NIR
# b5 <- raster('./NDVI/LC08_L1TP_036038_20130627_20180131_01_T1_B5.tif')
# # SWIR1
# b6 <- raster('./NDVI/LC08_L1TP_036038_20130627_20180131_01_T1_B6.tif')
# # SWIR2
# b7 <- raster('./NDVI/LC08_L1TP_036038_20130627_20180131_01_T1_B7.tif')
# # Panchromatic
# b8 <- raster('./NDVI/LC08_L1TP_036038_20130627_20180131_01_T1_B8.tif')
# # Cirrus
# b9 <- raster('./NDVI/LC08_L1TP_036038_20130627_20180131_01_T1_B9.tif')
# # Thermal Infrared (TIRS)1
# b10 <- raster('./NDVI/LC08_L1TP_036038_20130627_20180131_01_T1_B10.tif')
# # Thermal Infrared (TIRS)2
# b11 <- raster('./NDVI/LC08_L1TP_036038_20130627_20180131_01_T1_B11.tif')

# landsat.NDVI <- stack(b7, b6, b5, b4, b3, b2)
# landsat.NDVI

# Used the alredy cropped  images from before
# landsatRGB <- landsatcropA[[c(4,2,3)]]
# landsatFCC <- landsatcropB[[c(5,4,3)]]

###########################
## BASIC FUNCTION FOR NDVI:
library(rasterVis)

# img is a mutilayer Raster* object and i and k are the indices of the layers (layer numbers) 
# used to compute the vegetation index.
vi <- function(img, k, i) { 
  bk <- img[[k]]
  bi <- img[[i]]
  vi <- (bk - bi) / (bk + bi) 
  return(vi)
}

## For Landsat NIR = 5, red = 4. (bands 1,2 in object)
# STONE CANYON:

ndvi.SC <- vi(landsatcrop.B, 1, 2)
# plot(ndvi.SC, col = rev(terrain.colors(6)))
ndvi.raster.SC<-gplot(ndvi.SC) + geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradient(low = 'white', high = 'green') +
  coord_equal() +
  theme_bw() +
  labs(x="Easting (m)", y="Northing (m)") +
  theme(legend.position = c(.11,.73), legend.background = element_rect(colour = "black"),
        strip.text = element_blank())

# OWL HEAD:

ndvi.OH <- vi(landsatcrop.B2, 1, 2)
# plot(ndvi.OH, col = rev(terrain.colors(5)))

ndvi.raster.OH<-gplot(ndvi.OH) + geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradient(low = 'white', high = 'yellow') +
  coord_equal() +
  theme_bw() +
  labs(x="Easting (m)", y="Northing (m)") +
  theme(legend.position = c(.11,.73), legend.background = element_rect(colour = "black"),
        strip.text = element_blank())

ggarrange(ndvi.raster.OH, ndvi.raster.SC,  labels = c("A", "B"),
          ncol = 1, nrow = 2)

## Alternative function (didnt use this one)
# vi2 <- function(x, y) { 
#   (x - y) / (x + y)
# }

# ndvi2 <- overlay(landsat[[5]], landsat[[4]], fun=vi2) 
# plot(ndvi2, col=rev(terrain.colors(10)), main="Landsat-NDVI")

## We can explore the distribution of values contained within our raster using the 
## hist() function which produces a histogram. Histograms are often useful in identifying 
## outliers and bad data values in our raster data.

par(mfrow = c(2,2))

# view histogram of data
# STONE CANYON
Hist1<-hist(ndvi.SC,
     main = "",
     xlab = "NDVI",
     ylab= "Frequency",
     col = "wheat",
     xlim = c(0.05, 0.7),
     breaks = 30,
     xaxt = 'n') 
axis(side=1, at = seq(0.05,0.7, 0.05), labels = seq(0.05,0.7, 0.05))

# OWL HEAD
non.NDVI.hist <- hist(ndvi.OH,
     main = "",
     xlab = "NDVI",
     ylab= "Frequency",
     col = "wheat",
     xlim = c(0.05, 0.25),
     breaks = 30,
     xaxt = 'n')
non.NDVI.hist <- axis(side=1, at = seq(0.05, 0.25, 0.05), labels = seq(0.05, 0.25, 0.05))

# xlim = c(0.05, 0.7),
# breaks = 30,
# xaxt = 'n')
# axis(side=1, at = seq(0.05, 0.7, 0.05), labels = seq(0.05, 0.7, 0.05))




##############################################
##               THRESHOLDING
##############################################

## Can apply basic rules to get an estimate of spatial extent of different Earth surface 
## features. *Note: NDVI values are standardized and ranges between -1 to 1 (Higher 
## values indicate more green cover)

## Cells with NDVI values greater than 0.4 are definitely vegetation. Following operation 
## masks all cells that MIGHT NOT vegetation
# STONE CANYON
veg <- reclassify(ndvi.SC, cbind(-Inf, 0.3, NA)) 
plot(veg, main='Vegetation for SC')
# OWL HEADS
veg2 <- reclassify(ndvi.OH, cbind(-Inf, 0.3, NA)) 
plot(veg2, main='Vegetation for OH')

## Map area that corresponds to the peak between 0.25 and 0.3 in the NDVI histogram
# STONE CANYON
land <- reclassify(ndvi.SC, c(-Inf, 0.25, NA, 0.25, 0.3, 1, 0.3, Inf, NA)) 
plot(land, main = 'What is it? SC')
# OWL HEAD
land2 <- reclassify(ndvi.OH, c(-Inf, 0.25, NA, 0.25, 0.3, 1, 0.3, Inf, NA)) 
plot(land2, main = 'What is it? OH')

## Can plot land over original landsatFCC raster to find out more
plotRGB(landsatcropB, r=1, g=2, b=3, axes=TRUE, stretch="lin", 
        main="Landsat False Color Composite")
plot(veg, add=TRUE, legend=FALSE)

?`reclassify,Raster-method`

## Can create classes for different amounts of vegetation presence
#STONE CANYON
vegc <- reclassify(ndvi.SC, c(-Inf,0.25,1, 0.25,0.3,2, 0.3,0.4,3, 0.4,0.5,4, 0.5,Inf, 5)) 
plot(vegc,col = rev(terrain.colors(5)), main = 'NDVI based thresholding SC')
# OWL HEAD
vegc.OH <- reclassify(ndvi.OH, c(-Inf,0.25,1, 0.25,0.3,2, 0.3,0.4,3, 0.4,0.5,4, 0.5,Inf, 5)) 
plot(vegc.OH,col = rev(terrain.colors(5)), main = 'NDVI based thresholding OH')

# different classes perhaps more representative of Stone Canyon
# STONE CANYON
vegc.SC <- reclassify(ndvi.SC, c(-Inf,0.2,1, 0.2,0.25,2, 0.25,0.3,3, 0.3,0.4,4, 0.4,Inf, 5)) 
plot(vegc.SC,col = rev(terrain.colors(5)), main = 'NDVI based thresholding SC')
# OWL HEAD
vegc.OH <- reclassify(ndvi.OH, c(-Inf,-0.5,1, 0.09,0.1,2, 0.11,0.13,3, 0.3,0.6,4, 0.6,Inf, 5)) 
vegc.OH <- reclassify(ndvi.OH, c(-Inf,0.2,1, 0.2,0.25,2, 0.25,0.3,3, 0.3,0.4,4, 0.4,Inf, 5)) 
plot(vegc.OH,col = rev(terrain.colors(5)), main = 'NDVI based thresholding OH')




##############################################
##                  PCA
##############################################

## Multi-spectral data are sometimes transformed to helps to reduce the dimensionality and 
## noise in the data. The principal components transform is a generic data reduction method 
## that can be used to create a few uncorrelated bands from a larger set of correlated bands

## You can calculate the same number of principal components as the number of input bands. 
## The first principal component (PC) explains the largest percentage of variance and other
## PCs explain additional the variance in decreasing order

set.seed(1)
sr <- sampleRandom(landsat, 10000) plot(sr[,c(4,5)], main = "NIR-Red plot")

## This is known as vegetation and soil-line plot (Same as the scatter plot in earlier section)
pca <- prcomp(sr, scale = TRUE)
pca

pci <- predict(landsat, pca, index = 1:2)
plot(pci[[1]])

## The first principal component highlights the boundaries between land use classes or 
## spatial details, which is the most common information among all wavelengths. 
## Difficult to understand what the second principal component is highlighting. Try 
## thresholding again:

pc2 <- reclassify(pci[[2]], c(-Inf,0,1,0,Inf,NA))
par(mfrow = c(1,2))
plotRGB(landsatFCC, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin", 
        main = "Landsat False Color Composite")    
plotRGB(landsatFCC, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin", 
        main ="Landsat False Color Composite") plot(pc2, legend = FALSE, add = TRUE)

## ..... to be continued ..............
