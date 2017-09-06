setwd("/Users/rain/Suured_tegemised/Galmodel_phot/Input_sdss")
require(FITSio)
im<-readFITS(file=paste0("z.fits"))$imDat

box<-function(x0, y0, dx, dy, poleoluline){
  mask<-im==im
  x<-seq(round(x0-dx/2), round(x0+dx/2))
  y<-seq(round(y0-dy/2), round(y0+dy/2))
  x<-x[x>0 & x<dim(im)[1]]
  y<-y[y>0 & y<dim(im)[2]]
  mask[x,y]<-FALSE
  mask
}

mask<-box(51.5,324,13,22,0) * box(232.5,254.5,31,21,0) * box(290,270.5,16,11,0) * box(254.5,221.5,19,13,0) * box(303,84.5,36,37,0) * box(409.5,134.5,17,19,0) * box(228,93.5,10,9,0) * box(171,82,18,16,0) * box(231.5,7.5,25,17,0) * box(14.5,244,19,14,0) * box(379,140.5,14,9,0) * box(203.5,20.5,19,13,0)
writeFITSim(X = mask, file = "mask.fits")

mask[]<-0
writeFITSim(X = mask, file = "sigma.fits")