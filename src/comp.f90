module comp_module
! 	use constants
	use profile_collector_module
	use yldine_matemaatika_module
	
	type comp_type
		!kasutaja poolt sisestatud asjad
		character(len=default_character_length) 						:: comp_prof_name !profiili nimi ... nt "Einasto" vms
		character(len=default_character_length) 						:: comp_name !komponendi enda nimi ... lihtsustab testimisi jm
		character(len=default_character_length) 						:: comp_type_name !tyybi nimi ... valikud: stellar, attenuation, gas, DM
		character(len=default_character_length) 						:: population_name !tyybi nimi ... valikud: stellar, attenuation, gas, DM
		class(prof_den_base_type), allocatable	:: prof_den !profiil, mis saadakse profiles_den moodulist.... peab olema extended mingist pohityybist
		real(rk)			 					:: incl    !kaldenurk 0=faceon pi/2 = edge on
		real(rk)			 					:: dist    !galaktika kaugus kiloparsekites
		real(rk)								:: cnt_x   !komponendi tsenter vorreldes ruumi koordinaatidega
		real(rk)								:: cnt_y   !komponendi tsenter vorreldes ruumi koordinaatidega
		real(rk)								:: pos     !komponendi pos nurk vorreldes X teljega
		real(rk)								:: theta0  !komponendi k22ne taustsysteemi suhtes
		logical									:: pos_pool !milline pool on meie pool... defineerime kui pa+90 deg on meie pool
		!automaatselt m22ratavad
		real(rk) 								:: mass_abs_tol !see kasulik mudelpildi arvutamiseks, et liiga t2pselt ei teeks seda
		real(rk) 								:: sin_incl, cos_incl, tan_incl, sec_incl 
		real(rk) 								:: sin_pos, cos_pos, tan_pos, sec_pos
		integer 								:: comp_image_number !-1 on default negatiivne, et teeks uue pildi
	contains
		procedure, pass :: XpYp_to_XcYc => convert_XpYp_to_XcYc
	end type comp_type

	type prof_par_list_type !niisugust on vaja kuna ei tea mitu parameetrit vaja sisse lugeda
		logical :: filled = .false.
		type(par_type_real) :: par
		character(len=default_character_length) :: par_name
		type(prof_par_list_type), pointer :: next => null()
	end type prof_par_list_type
	!comp_input_type on sisendi sisselugemiseks... sisuliselt reaalarvud on prior piirkonnad ning fittimise linnikesed igal pool
	type comp_input_type
		character(len=default_character_length) 						:: comp_prof_name !profiili nimi ... nt "Einasto" vms
		character(len=default_character_length) 						:: comp_name !komponendi enda nimi ... lihtsustab testimisi jm
		character(len=default_character_length) 						:: comp_type_name !tyybi nimi ... valikud: stellar, attenuation, gas, DM
		character(len=default_character_length) 						:: population_name !tyybi nimi ... valikud: stellar, attenuation, gas, DM
		type(par_type_real)	 					:: incl !kaldenurk 0=faceon pi/2 = edge on
		type(par_type_real)			 					:: dist !galaktika kaugus kiloparsekites
		type(par_type_real)								:: cnt_x !komponendi tsenter vorreldes ruumi koordinaatidega
		type(par_type_real)								:: cnt_y !komponendi tsenter vorreldes ruumi koordinaatidega
		type(par_type_real)								:: pos !komponendi pos nurk vorreldes X teljega
		type(par_type_real)								:: theta0 !komponendi k22ne taustsysteemi suhtes
		type(prof_par_list_type)					:: prof_pars
		logical									:: pos_pool !milline pool on meie pool... defineerime kui pa+90 deg on meie pool
	contains
	end type comp_input_type
	
	contains
	
	subroutine init_comp(comp)
		implicit none
		type(comp_type), intent(inout) :: comp
! 			print*, "init comp juures"
		comp%sin_incl = sin(comp%incl)
		comp%cos_incl = cos(comp%incl)
		comp%tan_incl = tan(comp%incl)
		comp%sec_incl = 1.0/comp%cos_incl
		comp%sin_pos = sin(comp%pos)
		comp%cos_pos = cos(comp%pos)
		comp%tan_pos = tan(comp%pos)
		comp%sec_pos = 1.0/comp%cos_pos
	end subroutine init_comp
	elemental subroutine convert_XpYp_to_XcYc(comp, Xp, Yp, Xc, Yc)
		implicit none
		real(rk), intent(in) :: Xp, Yp
		real(rk), intent(out) :: Xc, Yc
		class(comp_type), intent(in) :: comp
		real(rk) :: x,y
		x = tan( (Xp-comp%cnt_x) ) * comp%dist
		y = tan( (Yp-comp%cnt_y) ) * comp%dist
		call coordinate_rotation(x, y, -1*comp%sin_pos, -1*comp%cos_pos, Xc, Yc)
	end subroutine convert_XpYp_to_XcYc
	elemental subroutine XcYcl_to_Rztheta(Xc, Yc, l,  sin_incl, cos_incl, tan_incl, sec_incl, theta0, R, z, theta) !Kontrollimata!
		implicit none
		real(rk), intent(in) :: Xc, Yc, l
		real(rk), intent(out) :: R,z,theta
		real(rk), intent(in) :: sin_incl, cos_incl, tan_incl, sec_incl, theta0
		z = l * cos_incl
		R = sqrt(Xc*Xc +  (z*tan_incl + Yc*sec_incl)**2)
		theta = sin_incl*Xc/R + theta0 !vaja KONTROLLIDA!
	end subroutine XcYcl_to_Rztheta
	
	
end module comp_module