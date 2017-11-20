module comp_module
! 	use constants
	use profile_collector_module
	use yldine_matemaatika_module
	use populations_module
	
	type comp_type
		!kasutaja poolt sisestatud asjad
		character(len=default_character_length) 						:: comp_prof_name !profiili nimi ... nt "Einasto" vms
		character(len=default_character_length) 						:: comp_name !komponendi enda nimi ... lihtsustab testimisi jm
		character(len=default_character_length) 						:: comp_type_name !tyybi nimi ... valikud: stellar, attenuation, gas, DM
		character(len=default_character_length) 						:: population_name !
		integer 														:: population_number !populations moodulist vastav number
		class(prof_den_base_type), allocatable	:: prof_den !profiil, mis saadakse profiles_den moodulist.... peab olema extended mingist pohityybist
		real(rk)			 					:: incl    !kaldenurk 0=faceon pi/2 = edge on
		real(rk)			 					:: dist    !galaktika kaugus kiloparsekites
		real(rk)								:: cnt_x   !komponendi tsenter vorreldes ruumi koordinaatidega
		real(rk)								:: cnt_y   !komponendi tsenter vorreldes ruumi koordinaatidega
		real(rk)								:: pos     !komponendi pos nurk vorreldes X teljega
		real(rk)								:: theta0  !komponendi k22ne taustsysteemi suhtes
		logical									:: pos_pool !milline pool on meie pool... defineerime kui pa+90 deg on meie pool
		!automaatselt m22ratavad
		real(rk) 								:: mass_abs_tol !see kasulik mudelpildi arvutamiseks, et liiga t2pselt ei teeks seda... seotud piksli suurusega
		real(rk)								:: massi_abs_tol_los !sirge peal olles t2psus
		real(rk) 								:: sin_incl, cos_incl, tan_incl, sec_incl 
		real(rk) 								:: sin_pos, cos_pos, tan_pos, sec_pos
		integer 								:: adaptive_image_number, adaptive_im_enne_tasandit, adaptive_im_p2rast_tasandit !-1 on default negatiivne, et teeks uue pildi
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
			if(.true.) then
				if(comp%incl<1.0e-5 ) print*, "incl liiga v2ike, s2tib suuremaks", comp%incl
				if(comp%incl<1.0e-5 ) comp%incl = 1.0e-5
				if(comp%incl>89.99*pi/180.0) print*, "incl liiga suur... s2tib v2iksemaks", comp%incl
				if(comp%incl>89.99*pi/180.0) comp%incl = 89.99*pi/180.0
			end if
		
		comp%sin_incl = sin(comp%incl)
		comp%cos_incl = cos(comp%incl)
		comp%tan_incl = tan(comp%incl)
		comp%sec_incl = 1.0/comp%cos_incl
		comp%sin_pos = sin(comp%pos)
		comp%cos_pos = cos(comp%pos)
		comp%tan_pos = tan(comp%pos)
		comp%sec_pos = 1.0/comp%cos_pos
		comp%population_number = get_pop_number_from_name(comp%population_name)
	contains
		function get_pop_number_from_name(pop_name) result(i)
			integer :: i
			character(len=default_character_length), intent(in) :: pop_name
			do i=1,size(populations, 1)
				if(trim(populations(i)%name) == trim(pop_name)) then
					return
				end if
			end do
			print*,  "Err: cannot match the population name with existing populations", trim(pop_name)
			stop
		end function
	end subroutine init_comp
	
	function leia_vabade_parameetrite_arv(input_comps) result(res)
		implicit none
		integer :: i
		type(prof_par_list_type), pointer ::  par_list
		type(comp_input_type), dimension(:), allocatable, target :: input_comps
		integer :: res
		res = 0
		do i=1,size(input_comps)
			if(input_comps(i)%incl%kas_fitib) res=res+1
			if(input_comps(i)%cnt_x%kas_fitib) res=res+1
			if(input_comps(i)%cnt_y%kas_fitib) res=res+1
			if(input_comps(i)%pos%kas_fitib) res=res+1
			if(input_comps(i)%theta0%kas_fitib) res=res+1
			par_list=>input_comps(i)%prof_pars
			do while(par_list%filled)
				if(par_list%par%kas_fitib) res=res+1
				if(associated(par_list%next)) then
					par_list => par_list%next
				else
					exit
				end if
			end do
		end do
		nullify(par_list)
	end function leia_vabade_parameetrite_arv
	
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
		R = sqrt(Xc*Xc +  (z*tan_incl - Yc*sec_incl)**2)
		theta = sin_incl*Xc/R + theta0 !TODO vaja KONTROLLIDA!
	end subroutine XcYcl_to_Rztheta
	
	
end module comp_module