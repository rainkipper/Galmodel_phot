module constants_module
	implicit none
	!progemise parameetrid
	integer, parameter :: rk = kind(1.0d0)
	integer, parameter :: default_character_length = 256
	integer, parameter :: max_prof_par_size = 1000 !ehk kokku fittimise peale tohib olla niipalju parameetreid
	!yldised konstandid
	real(rk), parameter :: pi = 4.0_rk*atan(1.0_rk)
	real(rk), parameter :: sqrt2 = sqrt(2.0_rk)
	!
	real(rk), parameter :: arcsec_to_rad = pi/180.0 / 60.0 / 60.0 !ehk sellega korrutades saab kaaresekundid radiaanidesse
	
	!muud parameetrid
	
	character(len=*), parameter :: model_filename_end = "_mdl.fits"
	character(len=*), parameter :: residual_filename_end = "_res.fits"
	character(len=*), parameter :: diff_filename_end = "_dif.fits"
	character(len=default_character_length)  :: multinest_output_header = "Output/MN_output/"
	character(len=default_character_length)  :: output_fit_file = "Output/out_parameters.txt"
	character(len=default_character_length)  :: output_ML_file = "Output/out_ML.txt"
	character(len=default_character_length)  :: sol_file_name
	character(len=default_character_length)  :: populations_file_name
	character(len=default_character_length)  :: filters_file_name
	logical :: kas_vaikselt = .false.
	
	!seaded
	integer  :: mis_fittimise_algoritm = 2
	integer  :: mis_fittimise_tyyp = 1 !mis tyypi fittimine... 1=populatsioonid ja massid,  2=komponendid ja ML
	logical  :: kas_fitib_massid_eraldi = .true. !iga likelihoodi arvutamise juures, kas fitib sisemiselt massid eraldi
	logical  :: kas_barrier = .true. !see, kas masside eraldi fittimisel on 0 juurde barj22r pandud
	logical  :: kas_koik_pildid_samast_vaatlusest = .true. 
	logical  :: via_adaptive_im = .true.
	logical  :: kas_los = .false. !kas profiili integreeritakse vaatejoont pidi... muidu votab ainult tasandis oleva v22rtuse 
	real(rk) :: max_calculation_time = 21600.0 !max arvutusaeg (s)... default on 6h
	logical  :: kas_rakendab_psf = .false. !
	logical  :: kas_psf_crop = .true. !kas kahandab psf suurust, et arvutused kiiremini oleks. 
	real(rk) :: psf_sisse_peab_j22ma = 0.99 !kui psf servasid natuke k2rbib, siis saab arvutamise kiiremaks
	integer  :: N_multinest_extra_points = 10
	real(rk) :: multinest_efr = 0.7
 	
	real(rk) :: default_los_kauguse_piir = 20.0_rk !integreerimise rada
	integer  :: pix_iter_maxlevel = 5
	real(rk) :: pix_edasi_jagamise_rel_t2psus = 0.005
	real(rk) :: abs_tol_kordaja = 0.01 !ehk see kordaja korrutatakse sigma-piltidest saadud kordajaga
	
	real(rk) :: massif_fiti_rel_t2psus = 0.003 !suhteline t2psus, mille korral loeb koondunuks masside eraldi fittimise
	real(rk) :: massi_fittimise_hyppe_kordaja = 0.7 !ehk kui kiiresti liigub iteratsioonide vahel
	real(rk) :: massi_fiti_lambda = 5.0 !m22rab kui t2pselt ei tohi massid nulli minna... pot barj22ris/chisq on see kui  -lambda*log(massi_kordaja)
	real(rk) :: lainepikkuse_dl = 0.1 !dl, yle mille integreeritakse
	
	real(rk) :: amoeba_fractional_tolerance = 0.01
	
	integer :: integration_min_iter=3, integration_max_iter=15 !integreerimise t2psuse parameetrid
	
	real(rk) :: adaptive_image_dist_piirang             = 1.0
	real(rk) :: adaptive_image_x0_default               = 0.0         !-35.0_rk
	real(rk) :: adaptive_image_y0_default               = 0.0         !-35.0_rk
	real(rk) :: adaptive_image_x1_default               =  35.0
	real(rk) :: adaptive_image_y1_default               =  35.0
	integer ::  adaptive_image_maxlevel                 = 12
	integer ::  adaptive_image_minlevel                 = 2
	real(rk) :: adaptive_image_edasijagamise_threshold  = 0.01            !suhteline jagamine
	real(rk) :: adaptive_image_min_spatial_resolution   = 0.005            !praegu ei lahuta rohkem kui 10pc

	logical :: kas_output_reana = .false.
	
	!tyybid jm
	type :: par_type !tyyp, mille eesm2rk on fititavaid parameetreid m2lus hoida. 
		integer  :: ref !millise komponendi v22rtust kasutab
		logical  :: kas_fitib !kas fititav parameeter
! 		logical  :: kas_par_muutus
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
	
contains
	subroutine read_constants_and_parameters(filename)
		implicit none
		character(len=default_character_length), intent(in) :: filename
		integer :: iunit, ios !molemad failindusega seotud muutujad
		integer :: id1, id2 !j2rjed sisseloetud sonedes
		character(len=default_character_length) :: rida, subrida, tyhi
		iunit = 22
		open(file = filename, unit = iunit, action = "read")
		ios = 1
		do while(ios .ge. 0)
			rida = tyhi !m2lu nullimine ymberkasutamiseks... selleks, et trim k2suga asjade allaj22nud asju ei kasutaks
			read(fmt = "(a)", iostat=ios, unit=iunit) rida
			rida = adjustl(rida) !eemdalab koik esimesed tyhikud ja paneb loppu... 
			if( len(trim(rida))<2 .or. rida(1:1)=="#") cycle !mottetute ridade v2ltimine
			id1 = 1; id2 = index(rida, "=")
			if(id2==0) cycle
			subrida = tyhi
			subrida = trim(ADJUSTL(rida((id2+1):len(rida))))
			select case(trim(adjustl(rida(1:(id2-1)))))
			case("kas_fitib_massid_eraldi") ; read(subrida, fmt=*) kas_fitib_massid_eraldi
			case("kas_barrier"); read(subrida, fmt = *) kas_barrier
			case("kas_koik_pildid_samast_vaatlusest"); read(subrida, fmt = *) kas_koik_pildid_samast_vaatlusest
			case("via_adaptive_im"); read(subrida, fmt = *) via_adaptive_im
			case("kas_los"); read(subrida, fmt = *) kas_los
			case("kas_rakendab_psf"); read(subrida, fmt = *) kas_rakendab_psf
			case("kas_psf_crop"); read(subrida, fmt = *) kas_psf_crop
			case("psf_sisse_peab_j22ma"); read(subrida, fmt = *) psf_sisse_peab_j22ma
			case("pix_iter_maxlevel"); read(subrida, fmt = *) pix_iter_maxlevel
			case("massif_fiti_rel_t2psus"); read(subrida, fmt = *) massif_fiti_rel_t2psus
			case("pix_edasi_jagamise_rel_t2psus"); read(subrida, fmt = *) pix_edasi_jagamise_rel_t2psus	
			case("adaptive_image_x0_default"); read(subrida, fmt=*) adaptive_image_x0_default               
			case("adaptive_image_y0_default"); read(subrida, fmt=*) adaptive_image_y0_default               
			case("adaptive_image_x1_default"); read(subrida, fmt=*) adaptive_image_x1_default               
			case("adaptive_image_y1_default"); read(subrida, fmt=*) adaptive_image_y1_default               
			case("adaptive_image_maxlevel"); read(subrida, fmt=*) adaptive_image_maxlevel                 
			case("adaptive_image_minlevel"); read(subrida, fmt=*) adaptive_image_minlevel                 
			case("adaptive_image_edasijagamise_threshold"); read(subrida, fmt=*) adaptive_image_edasijagamise_threshold  
			case("adaptive_image_min_spatial_resolution"); read(subrida, fmt=*) adaptive_image_min_spatial_resolution   
			case("mis_fittimise_algoritm"); read(subrida, fmt=*)  mis_fittimise_algoritm
			case("massi_fiti_lambda"); read(subrida, fmt=*) massi_fiti_lambda
			case("massi_fittimise_hyppe_kordaja"); read(subrida, fmt=*) massi_fittimise_hyppe_kordaja
			case("amoeba_fractional_tolerance"); read(subrida, fmt=*) amoeba_fractional_tolerance
			case("mis_fittimise_tyyp"); read(subrida, fmt=*) mis_fittimise_tyyp
			case("default_los_kauguse_piir"); read(subrida, fmt=*) default_los_kauguse_piir
			case("abs_tol_kordaja"); read(subrida, fmt=*) abs_tol_kordaja
			case("adaptive_image_dist_piirang"); read(subrida, fmt=*) adaptive_image_dist_piirang
			case("integration_min_iter"); read(subrida, fmt=*) integration_min_iter
			case("integration_max_iter"); read(subrida, fmt=*) integration_max_iter
			case("N_multinest_extra_points"); read(subrida, fmt=*) N_multinest_extra_points	
			case("multinest_efr"); read(subrida, fmt=*) multinest_efr	
			case("multinest_output_header"); read(subrida, fmt="(a)") multinest_output_header	
			case("output_ML_file"); read(subrida, fmt="(a)") output_ML_file	
			case("output_fit_file"); read(subrida, fmt="(a)") output_fit_file	
			case("kas_vaikselt"); read(subrida, fmt = *) kas_vaikselt	
			case("kas_output_reana"); read(subrida, fmt = *) kas_output_reana	
			case("lainepikkuse_dl"); read(subrida, fmt = *) lainepikkuse_dl
			case("sol_file"); read(subrida, fmt = "(a)") sol_file_name
			case("populations_file"); read(subrida, fmt = "(a)") populations_file_name
			case("filters_file"); read(subrida, fmt = "(a)") filters_file_name
			case("max_calculation_time"); read(subrida, fmt = *) max_calculation_time
		
		
		
			
			case default
				print*, "Tundmatu parameeter sisselugemisel: ",trim(adjustl(rida(1:(id2-1))))
				stop
			end select
		end do
		close(unit = iunit)
		
	end subroutine read_constants_and_parameters
	
		
	
end module constants_module