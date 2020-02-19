source("/home/wietse/Documents/RProjects/metGeneratoR/melo_functions.R")
###########################
dates <- seq(as.POSIXct('2014-01-01 00:00:00'), by = "hours", length = 24)

# potential_radiation(date=c('2014-01-01 12:00:00','2014-01-01 13:00:00'),lon=8.86,lat= 51.00,timezone=1)
# rad_hour <- disaggregate_radiation(radiation=10.6625,date=dates,lon=8.86,lat= 51.00,timezone=1)

### Read netcdf
library(ncdf4)
ncid<-nc_open("~/sw_1day.nc")
rsds<-ncvar_get(ncid,"SWdown")
nc_close(ncid)

###
radlats<-array(NA, dim=c(720, 360, 24))
iy<-1
for (ilat in seq(-90,89.75,.5)) {
  ix<-1
  print(iy)
  for (ilon in seq(-180,179.75,.5)) {
    radlats[ix,iy,] <- disaggregate_radiation(radiation=rsds[ix,iy],date=dates,lon=ilon,lat= ilat,timezone=2)
    # radlats[ix,iy,] <- potential_radiation(dates=dates,lon=ilon,lat= ilat,timezone=1)
    ix<-ix+1
  }
  iy<-iy+1
}

ncid<-nc_open("~/sw_3h.nc")
rsds_3h<-ncvar_get(ncid,"SWdown")
nc_close(ncid)
image(rsds_3h[,,1])

mask <- rsds_3h[,,1]
mask[mask>=0] <- 1
#image(mask)
max(rsds_3h[,,1], na.rm = T)
max(radlats[,,1], na.rm = T)


inData <- radlats
nhourly<-3
outData <- array(data = 0, dim = c(720,360,(24/nhourly)))
inrecs <- c(1:24)
outrecs <- rep(1:(24/nhourly), each = nhourly)

for(i in 1:24) outData[, , outrecs[i]] <- outData[, , outrecs[i]] + inData[, , inrecs[i]]
for(i in 1:8) outData[, , i] <- (outData[, , i] / 3 ) * mask


qq<-seq(1,24,3)
outData2<-outData
outData3<-outData
outData4<-outData
for (i in 1:8) {
  outData2[,,i] <- radlats[,,qq[i]]
  outData3[,,i] <- radlats[,,(qq[i]+1)]
  outData4[,,i] <- radlats[,,(qq[i]+2)]
}

ilat <- 200
ilon <- 360
pdf(paste0("rplot_",2,".pdf") )
par(mfrow=c(2,2))
plot(rsds_3h[ilon,ilat,])
lines(outData[ilon,ilat,])
lines(outData2[ilon,ilat,], col = "green")
lines(outData3[ilon,ilat,], col = "red")
lines(outData4[ilon,ilat,], col = "blue")
image(rsds_3h[,,4])
image(outData[,,4])
image(radlats[,,11]*mask)
max(rsds_3h[,,4], na.rm = T)
max(outData[,,4], na.rm = T)
max(radlats[,,11]*mask, na.rm = T)
# Close the pdf file
dev.off() 

