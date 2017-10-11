setwd("~/Github/Galmodel_phot/")
incl<-75.2 * pi/180 #kaldenurk
galdist<-10000 #kpc
lopmatus<-50
#pos nurka ei arvesta


Q_fun<-function(q,incl){
   sqrt( cos(incl)^2 + q^2*sin(incl)^2 )
}
surf_lum<-function(incl, func, A){
   e<-environment(func)
   integreeritav<-function(a) func(a)*a/sqrt(a^2-A^2)
   tmp<-integrate(integreeritav, lower = A, upper = ifelse(A>(lopmatus/2), 3*A, lopmatus))$value
   Q<-Q_fun(e$q, incl)
   2*e$q/Q*tmp
}
gen_Einasto1<-function(ac, q, N, dN, rhoc){
   force(N); force(q); force(ac); force(dN); force(rhoc)
   function(a){
      #a<-sqrt(R^2+(z/q)^2)
      rhoc*exp(-dN*((a/ac)^(1/N) - 1))
   }
}

yhiku_kordaja<-0.1
Nucleus<-gen_Einasto1(ac = 0.0234, q = 0.99, N = 4.0, rhoc = yhiku_kordaja*1.713, dN = 11.668) 
Bulge<-gen_Einasto1(ac = 1.155, q = 0.72, N = 2.7, dN = 7.769, rhoc = yhiku_kordaja*9.201e-1)
Disc<-gen_Einasto1(ac = 10.67, q = 0.17, N = 1.2, dN = 3.273, rhoc = yhiku_kordaja*1.307e-2)
Young_disc<-gen_Einasto1(ac = 11.83, q = 0.01, N = 0.2, dN = 0.316, rhoc = yhiku_kordaja*1.179e-2) 
Stellar_halo<-gen_Einasto1(ac = 12.22, q = 0.50, N = 3.0, dN = 8.669, rhoc = yhiku_kordaja*4.459e-4)
mask<-2:5
funcs<-list(Nucleus, Bulge, Disc, Young_disc, Stellar_halo)[mask]
Q_list<-unlist(Map(Q_fun, c(0.99, 0.72,0.17,0.01,0.50)[mask], incl ))
ML_g<-c(4.44, 5.34, 5.23, 1.23, 6.19)[mask]
ML_r<-c(3.2, 4.08, 3.92, 1.12, 4.48)[mask]
ML_i<-c(2.35, 3.01, 2.92, 0.88, 3.25)[mask]

pildi_piir<-20; N<-100
coord<-seq(-pildi_piir, pildi_piir, length.out = N)
grid<-expand.grid(X = coord, Y = coord);
grid$lum<-0

for(i in seq_along(funcs)){
   grid$A<-sqrt(grid$X^2 + (grid$Y/Q_list[i])^2)
   tmp<-Vectorize(surf_lum, vectorize.args = "A")
   grid[[paste0("C",i)]]<- tmp(incl = incl, func = funcs[[i]], A = grid$A)
}



filters<-c("g", "r", "i")

for(f in seq_along(filters)){
   L<-grid$X*0 #algne
   for(i in seq_along(ML_g)) L <- L + get(paste0("ML_",filters[f]))[i]*grid[[paste0("C",i)]]
   writeFITSim(X = matrix(L, ncol=N, nrow=N), file = paste0("Input/Mock/",filters[f],".fits"))
}
writeFITSim(X = matrix(L*0+1, ncol=N, nrow=N), file = paste0("Input/Mock/mask.fits"))

#psf tegemine
N_psf = 7
psf<-matrix(0, N_psf,N_psf); cnt_val<-0.8
tmp<-(N_psf-1)/2+1
psf[]<-(1-cnt_val)/(length(psf)-1); psf[tmp,tmp]<-cnt_val
writeFITSim(psf, file = paste0("Input/Mock/psf.fits"))

scale<-(coord[2]-coord[1])/galdist*180/pi*3600
image(log(mx), useRaster = TRUE, asp=TRUE, x=coord, y = coord)


print(paste("Dist =", galdist, "kpc"))
print(paste("1 pixel is ", scale ," arcsec"))















