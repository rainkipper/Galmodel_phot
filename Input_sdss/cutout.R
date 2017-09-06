setwd("/Users/rain/Suured_tegemised/Galmodel_phot/Input_sdss")
require(FITSio)
filtrid<-c("u", "g", "r", "i", "z")

par_to_mask<-function(im, x0, dx, y0, dy){
  mask<-im!=im
  x<-seq(round(x0-dx/2), round(x0+dx/2))
  y<-seq(round(y0-dy/2), round(y0+dy/2))
  x<-x[x>0 & x<dim(im)[1]]
  y<-y[y>0 & y<dim(im)[2]]
  mask[x,y]<-TRUE
  mask
}

for(filter in filtrid){
  im<-readFITS(file=paste0("frame-",filter,"-004135-6-0288.fits"))$imDat
  mask<-par_to_mask(im, 803, 412, 794, 347)
  x<-seq(round(803-412/2), round(803+412/2))
  y<-seq(round(794-347/2), round(794+347/2))
  writeFITSim(X = im[x,y], file=paste0(filter, ".fits"))
}
