#r filter baasiks... fitib k2sitsi
setwd("/Users/rain/Suured_tegemised/Galmodel_phot/Input_sdss")
require(FITSio)
if(FALSE){
  u<-readFITS("u.fits")$imDat; 
  g<-readFITS("g.fits")$imDat; 
  r<-readFITS("r.fits")$imDat; 
  i<-readFITS("i.fits")$imDat; 
  z<-readFITS("z.fits")$imDat; 
}
nihe<-function(mx, dx, dy){
  res<-mx-mx
  N<-dim(mx)[1]
  for(i in 1:N){
    j<-i+dx
    if(j<0.5) j<-round(abs(j+N))
    if(j>(N+0.5)) j<-round(j-N)
    res[i,]<-mx[j,]
  }
  mx<-res
  N<-dim(mx)[2]
    for(i in 1:N){
      j<-i+dy
      if(j<0.5) j<-round(abs(j+N))
      if(j>(N+0.5)) j<-round(j-N)
      res[,i]<-mx[,j]
    }
  res
}
pilt<-function(mx) image(log(mx - 2*min(mx)), useRaster = TRUE)

# mx<-nihe(z, 2, 8)
# mx<-nihe(i, 2, 2)
# mx<-nihe(g, 3, 12)
# mx<-nihe(u, 5, 5)
# pilt(r); lines(contour(log(mx+2*min(mx)), nlevels = 5, add=TRUE), col="blue")
if(FALSE){
  writeFITSim(nihe(u, 5, 5), "u_nihutatud.fits")
  writeFITSim(nihe(g, 3, 12), "g_nihutatud.fits")
  writeFITSim(nihe(i, 2, 2), "i_nihutatud.fits")
  writeFITSim(nihe(z, 2, 8), "z_nihutatud.fits")
}