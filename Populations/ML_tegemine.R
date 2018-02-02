setwd("/Users/rain/Suured_tegemised/Galmodel_phot/Populations")
# setwd("/Users/rain/Suured_tegemised/Califa_SDSS/Populatsioonide_vidin")
filters<-c("u", "g", "r", "i", "z")

if(TRUE){
  t<-read.table("blanton_spec_raw.txt", as.is=TRUE)
  #For each of the five templates, we have the template spectrum
  # spec in units of ergs s-1 cm-2 A-1 Msun-1, as it would be ob-  
  # served at 10 pc distance.
  names(t)<-c("lambda", "T1", "T2", "T3", "T4", "T5")
  
  sol<-read.table("solar_standard_spectrum.csv", header = FALSE, as.is=TRUE, comment.char = "#")
  #W*m-2*nm-1
  names(sol)<-c("lambda", "lum1", "tilt", "lum2")
  sol_unit_to_spec_unit<-1e2 #pohiyhikud erg yhikuteks
  sol_unit_to_spec_unit<-sol_unit_to_spec_unit * tan(1.0/3600*pi/180/10)^2 #kaugus teisendub
  sol$lambda<-sol$lambda*10 #koik angstromideks
  sol$lum1<-sol$lum1 * sol_unit_to_spec_unit
  sol$lum2<-sol$lum2 * sol_unit_to_spec_unit
} #andmete lugemine ja filtri lugemine
sol_fun<-approxfun(sol$lum1~sol$lambda, rule = 2)
cat("sol"); print(range(sol$lambda))
cat("spec"); print(range(t$lambda))

res<-numeric(length(filters)*5)*NA
res<-data.frame(filt=res, template=res, lum_sum=res, ML=res)
counter<-0
kas_tolmuga<-FALSE


for(filt_name in filters){   
  f<-read.table(paste0("",filt_name,".dat.txt"), as.is=TRUE)
  names(f)<-c("lambda", "T_point", "T_airmass", "T", "T_extinction")
  filter_fun<-approxfun(f$T_point~f$lambda)
  cat("filt"); print(range(f$lambda))
  for(spec_nr in 1:5){
    counter<-counter+1
    spec_fun<-approxfun(t[[paste0("T",spec_nr)]]~t$lambda)
    
    if(TRUE){
      dl<-0.1
      l<-seq(min(f$lambda),max(f$lambda),dl)
      specsum<-sum(spec_fun(l)*filter_fun(l))*dl
      solsum<-sum(sol_fun(l)*filter_fun(l))*dl
      # print(paste(filt_name, spec_nr, specsum/solsum, solsum/specsum))
      
      res$filt[counter]<-filt_name
      res$template[counter]<-spec_nr
      res$lum_sum[counter]<-specsum/solsum
      res$ML[counter]<-solsum/specsum
    }#ML arvutamine
    
  }}
if(kas_tolmuga){
}else{
  write.table(res, file = "ML.txt", row.names=FALSE, quote = FALSE)
}


if(!TRUE){
  tekita_filt<-function(filt_name){
    f<-read.table(paste0("",filt_name,".dat.txt"), as.is=TRUE)
    names(f)<-c("lambda", "T_point", "T_airmass", "T", "T_extinction")
    approxfun(f$T_point~f$lambda)
  }
  flts<-Map(tekita_filt, filters)
  
  pdf("templates_and_fitlers.pdf", width=8.4/2.54*2, 8.4/2.54*2.4, pointsize = 7)
  par(mfrow=c(3,2))
  require(magicaxis)
  v2rvid<-c("blue", "darkgreen", "green", "darkred", "red")
  mask<-t$lambda>3000 & t$lambda<9000
  plot(0,0, xlim=c(2000, 9000), ylim=c(0,3e-3), xlab="lambda (A)", ylab="", xaxt="n", yaxt="n", main="Just templates")
  f<-sum
  with(t[mask,], {
    lines(lambda, T1/f(T1), col=v2rvid[1])
    lines(lambda, T2/f(T2), col=v2rvid[2])
    lines(lambda, T3/f(T3), col=v2rvid[3])
    lines(lambda, T4/f(T4), col=v2rvid[4])
    lines(lambda, T5/f(T5), col=v2rvid[5])
  })
  magaxis(1:4, labels = c(TRUE,FALSE,FALSE,FALSE))
  legend("topleft", col=v2rvid, legend=paste("Blanton",1:5), pch=16, bg="transparent", box.lty=0)
  
  for(i in 1:length(filters)){
  plot(0,0, xlim=c(2000, 9000), ylim=c(0,3e-3), xlab="lambda (A)", ylab="", xaxt="n", yaxt="n", main=paste("SDSS",filters[i]))
  f<-sum
  with(t[mask,], {
    lines(lambda, flts[[i]](lambda)*T1/f(T1), col=v2rvid[1])
    lines(lambda, flts[[i]](lambda)*T2/f(T2), col=v2rvid[2])
    lines(lambda, flts[[i]](lambda)*T3/f(T3), col=v2rvid[3])
    lines(lambda, flts[[i]](lambda)*T4/f(T4), col=v2rvid[4])
    lines(lambda, flts[[i]](lambda)*T5/f(T5), col=v2rvid[5])
  })
  magaxis(1:4, labels = c(TRUE,FALSE,FALSE,FALSE))
  legend("topleft", col=v2rvid, legend=paste("Blanton",1:5), pch=16, bg="transparent", box.lty=0)
  }
  dev.off()
  
}#template joonised


