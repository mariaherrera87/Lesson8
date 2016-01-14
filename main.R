## Author : Maria Yi
## Jan 13,2016

library(sp)
library(raster)

load("data/GewataB1.rda")
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")
load("data/GewataB5.rda")
load("data/GewataB7.rda")
load("data/vcfGewata.rda")
## band 6 (thermal infra-red) will be excluded from this exercise
vcfGewata[vcfGewata > 100] <- NA
summary(vcfGewata)
alldata <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7, vcfGewata)
names(alldata) <- c("band1", "band2", "band3", "band4", "band5", "band7", "VCF")

## extract all data to a data.frame
df <- as.data.frame(getValues(alldata))
# define a model to process the VCF
lmodel <- lm(VCF ~ band7 + band5 + band4 + band3 + band2 + band1, data = df)
summary(lmodel)
dataprd <- predict(alldata, model = lmodel)

# plot the result of the two VCF
opar <- par(mfrow = c(1,2))
plot(vcfGewata,main = "Original VCF")
plot(dataprd, zlim=c(0,100), main="Predict VCF")

RMSE <- overlay(alldata$VCF, dataprd, fun = function(x,y){(x - y)^2})
RMSE <- sqrt(cellStats(RMSE, stat='mean'))

#---------------------RMSE for every training polygon------------------#
# we just use the predicted rasterlayer from RandomForest model,and save it as 'RFpred.grd'
library(rgdal)
getwd()
RFpred <- raster('data/RFpred.grd')
plot(RFpred,col = c("red","blue","green"))
pre <- freq(RFpred)
# get zonal stastistic
zonal(RFpred,z)

# Make an NA-value raster based on the LC raster attributes
cropmask <-formask <- wetmask <- setValues(raster(RFpred), NA)

# define the area for all corresponding classes
cropmask[RFpred==1] <- 1
formask[RFpred==2] <- 1
wetmask[RFpred==3] <- 1

# Mask original and predicted VCF raster and plot the values of different classes
par(mfrow = c(3,2))
cropdata <- mask(dataprd,cropmask)
croppred <- mask(alldata$VCF,cropmask)


fordata <- mask(dataprd,formask)
forpred <- mask(alldata$VCF,formask)

wetdata <- mask(dataprd,wetmask)
wetpred <- mask(alldata$VCF,wetmask)

#calculate RMSE for the three differnt classes
cropRMSE <- overlay(cropdata, croppred, fun = function(x,y){(x - y)^2})
cropRMSE <- sqrt(cellStats(cropRMSE, stat='mean', na.rm =T))

forRMSE <- overlay(fordata, forpred, fun = function(x,y){(x - y)^2})
forRMSE <- sqrt(cellStats(forRMSE, stat='mean', na.rm =T))

wetRMSE <- overlay(wetdata, wetpred, fun = function(x,y){(x - y)^2})
wetRMSE <- sqrt(cellStats(wetRMSE, stat='mean', na.rm =T))


# plot and compare difference between original and predict VCF
par(mfrow = c(3,2))
plot(cropdata,zlim = c(1,100),main = "Cropland Original VCF")
plot(croppred,zlim = c(1,100),main = "Cropland Predict VCF", xlab = paste("RMSE of crop is :",round(cropRMSE,digits = 1)))

plot(fordata,zlim = c(1,100),main = "Forest Original VCF")
plot(forpred,zlim = c(1,100),main = "Forest Predict VCF", xlab = paste("RMSE of forest is :",round(forRMSE,digits = 1)))

plot(wetdata,zlim = c(1,100),main = "Wetland Original VCF")
plot(wetpred,zlim = c(1,100),main = "Wetland Predict VCF", xlab = paste("RMSE of wetland is :",round(wetRMSE,digits = 1)))
