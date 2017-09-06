#http://www.bruzual.org/bc03/Updated_version_2016/ on saadud BC absmagid
Mag_sun<-c(6.45, 5.14, 4.65, 4.54, 4.52)
Mag_to_ML<-function(Mag){
  10**((Mag-Mag_sun)/5)
}
Mag_logT5<-c(2.4750,    2.4358 ,   2.8318 ,   3.1126  ,  3.3600)
Mag_logT6<-c(2.4589,    2.4232,    2.8203,    3.1019,    3.3500)
Mag_logT7<-c(2.0103,    2.0307 ,   2.4484 ,   2.7444   , 3.0057)
Mag_logT8<-c(4.2940,    3.5918,    3.7643 ,   3.6554   , 3.0769)
Mag_logT9<-c(7.3572,    5.9000,    5.3087 ,   4.9367  ,  4.4222)
Mag_logT10<-c(9.1697,    7.8893,    7.2510 ,   6.8825 ,   6.3934)

ML_6<-Mag_to_ML(Mag_logT6); print(ML_6)
ML_9<-Mag_to_ML(Mag_logT9); print(ML_9)
ML_10<-Mag_to_ML(Mag_logT10); print(ML_10)
