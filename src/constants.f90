module constants_module
	implicit none
	!progemise parameetrid
	integer, parameter :: rk = kind(1.0d0)
	integer, parameter :: default_character_length = 128
	integer, parameter :: max_prof_par_size = 1000 !ehk kokku fittimise peale tohib olla niipalju parameetreid
	!yldised konstandid
	real(rk), parameter :: pi = 4.0_rk*atan(1.0_rk)
	!
	real(rk), parameter :: arcsec_to_rad = pi/180.0 / 60.0 / 60.0 !ehk sellega korrutades saab kaaresekundid radiaanidesse
	
	character(len=*), parameter :: model_filename_end = "_mdl.fits"
	character(len=*), parameter :: residual_filename_end = "_res.fits"
	character(len=*), parameter :: diff_filename_end = "_dif.fits"
	
	!tyybid jm
	type :: par_type !tyyp, mille eesm2rk on fititavaid parameetreid m2lus hoida. 
		integer  :: ref !millise komponendi v22rtust kasutab
		logical  :: kas_fitib !kas fititav parameeter
		character(len=default_character_length), allocatable :: prior_type !priori tyyp
	end type par_type
	type, extends(par_type) :: par_type_real
		real(rk) :: val !reaalne v22rtus
		real(rk) :: min !priori miinimum
		real(rk) :: max !priori maksimum
	end type par_type_real
	type, extends(par_type) :: par_type_int
		integer  :: val !reaalne v22rtus
		integer  :: min !priori miinimum
		integer  :: max !priori maksimum
	end type par_type_int
	type, extends(par_type) :: par_type_logical
		logical  :: val !reaalne v22rtus
		logical  :: min !priori miinimum
		logical  :: max !priori maksimum
	end type par_type_logical
	
	
		
	
end module constants_module