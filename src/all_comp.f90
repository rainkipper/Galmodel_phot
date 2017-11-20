module all_comp_module
! 	use constants
	use comp_module

	
	type all_comp_type
		class(comp_type), allocatable, dimension(:) 		:: comp
		class(comp_type), allocatable, dimension(:)			:: dust_comp
		integer, allocatable, dimension(:) 					:: comp_im_ref !komponentide pildid... 2D massiiv kuna pilte voib olla mitut tyypi
		logical, allocatable, dimension(:)  				:: recalc_adaptive_ims
		integer 											:: N_comp !kogu komponentide arv		
		!accounting .. 											
		integer 											:: N_dust
	end type all_comp_type
	
! 	private type(all_comp_type) :: all_comp
	
contains
	!function fun_all_potential
	!function fun_all_den
	!function fun_dpot_dR
	!function fun_dpot_dz
	!function fun_dpot_dtheta

	subroutine convert_input_comp_to_all_comp(input_comps, all_comp)
		implicit none
		type(comp_input_type), intent(in), dimension(:), target :: input_comps
		type(all_comp_type), intent(inout) :: all_comp
		class(comp_type), pointer :: point_to_comp
		type(prof_par_list_type), pointer ::  par_list
		integer :: i
		
		call loenda_komponentide_tyypide_arvud(input_comps, all_comp%N_comp, all_comp%N_dust)
		if(not(allocated(all_comp%comp)) .and. all_comp%N_comp>0) allocate(all_comp%comp(1:all_comp%N_comp))
		if(not(allocated(all_comp%dust_comp)) .and. all_comp%N_dust>0) allocate(all_comp%dust_comp(1:all_comp%N_dust))
		do i=1,size(input_comps, 1)

			call point_to_right_component(i, point_to_comp, input_comps, all_comp)
			point_to_comp%comp_prof_name 	= input_comps(i)%comp_prof_name
			point_to_comp%comp_name 		= input_comps(i)%comp_name
			point_to_comp%comp_type_name  	= input_comps(i)%comp_type_name 
			point_to_comp%population_name 	= input_comps(i)%population_name
			point_to_comp%incl     			= input_comps(i)%incl%val    
			point_to_comp%dist     			= input_comps(i)%dist%val    
			point_to_comp%cnt_x    			= input_comps(i)%cnt_x%val   
			point_to_comp%cnt_y    			= input_comps(i)%cnt_y%val   
			point_to_comp%pos      			= input_comps(i)%pos%val     
			point_to_comp%theta0   			= input_comps(i)%theta0%val
			!
			!profiili allokeerimine ja parameetrite paika panemine... siia peab iga profiili korral tulema kontroll
			!
			if(trim(point_to_comp%comp_prof_name) == "Einasto" .and. not(allocated(point_to_comp%prof_den))) allocate(prof_Einasto_type::point_to_comp%prof_den)
			if(trim(point_to_comp%comp_prof_name) ==  "Sersic" .and. not(allocated(point_to_comp%prof_den))) allocate( prof_Sersic_type::point_to_comp%prof_den)
			if(trim(point_to_comp%comp_prof_name) ==  "dustplane" .and. not(allocated(point_to_comp%prof_den))) allocate( prof_dustplane_type::point_to_comp%prof_den)
			!
			! 
			!
			par_list=>input_comps(i)%prof_pars
			do while(par_list%filled)
				call point_to_comp%prof_den%set_val(trim(par_list%par_name), par_list%par%val)
				if(associated(par_list%next)) then
					par_list => par_list%next
				else
					exit
				end if
			end do
			nullify(par_list)
		end do
	end subroutine convert_input_comp_to_all_comp
	
	subroutine asenda_viited(input_comps, all_comp, recalc_comp)
		implicit none
		type(comp_input_type), intent(in), dimension(:), allocatable, target :: input_comps
		type(all_comp_type), intent(inout), target :: all_comp
		logical, intent(inout), dimension(:), allocatable, optional :: recalc_comp !siit lisatakse viidete kaudu ymberarvutuse vajadus
		type(prof_par_list_type), pointer ::  par_list
		integer :: i!,j
		integer :: viide
		real(rk) :: tmp!, tmp2
		logical :: kas_recalc !ehk kas on vaja t2ita ymber arvutuste asju
		integer, dimension(:), allocatable :: input_comps_ref_to_lum_comp_ref, input_comps_ref_to_dust_comp_ref
		class(comp_type), pointer :: comp_point, kuhu_viited_pannakse
			
		kas_recalc = present(recalc_comp) 
		
		do i=1,size(input_comps,1)
			call point_to_right_component(i, kuhu_viited_pannakse, input_comps, all_comp)
			call point_to_right_component(input_comps(i)%incl%ref, comp_point, input_comps, all_comp)
			kuhu_viited_pannakse%incl     = comp_point%incl
			call point_to_right_component(input_comps(i)%cnt_x%ref, comp_point, input_comps, all_comp)
			kuhu_viited_pannakse%cnt_x     = comp_point%cnt_x
			call point_to_right_component(input_comps(i)%cnt_y%ref, comp_point, input_comps, all_comp)
			kuhu_viited_pannakse%cnt_y     = comp_point%cnt_y
			call point_to_right_component(input_comps(i)%pos%ref, comp_point, input_comps, all_comp)
			kuhu_viited_pannakse%pos     = comp_point%pos
			call point_to_right_component(input_comps(i)%theta0%ref, comp_point, input_comps, all_comp)
			kuhu_viited_pannakse%theta0     = comp_point%theta0

			par_list=>input_comps(i)%prof_pars
			do while(par_list%filled)
				call point_to_right_component(par_list%par%ref, comp_point, input_comps, all_comp)
				call comp_point%prof_den%get_val(trim(par_list%par_name), tmp) !votab v22rtuse ref jaoks
				call kuhu_viited_pannakse%prof_den%set_val(trim(par_list%par_name), tmp) !paneb v22rtuse teise kohta paika
				if(associated(par_list%next)) then
					par_list => par_list%next
				else
					exit
				end if
			end do
			call kuhu_viited_pannakse%prof_den%init_profile()
			call init_comp(kuhu_viited_pannakse)
		end do
		nullify(par_list)
	contains
	end subroutine asenda_viited
	subroutine point_to_right_component(input_comp_number, res, input_comps, all_comp)
		implicit none
		integer, intent(in) :: input_comp_number
		type(comp_input_type), intent(in), dimension(:), target :: input_comps
		type(all_comp_type), intent(in), target :: all_comp
		class(comp_type), intent(out), pointer :: res
		integer :: i, Nlum, Ndust
		Nlum = 0; Ndust=0
		do i=1,input_comp_number
			select case (trim(input_comps(i)%comp_type_name))
			case("stellar"); Nlum = Nlum + 1
			case("dust"); Ndust = Ndust + 1
			case default
				print*, "unknown type of component in all_comp.f90:", trim(input_comps(i)%comp_type_name)
				stop 
			end select
		end do
		select case (trim(input_comps(input_comp_number)%comp_type_name))
		case("stellar"); res => all_comp%comp(Nlum)
		case("dust"); res => all_comp%dust_comp(Ndust)
		end select
	end subroutine point_to_right_component
	subroutine loenda_komponentide_tyypide_arvud(input_comps, Nlum, Ndust)
		implicit none
		type(comp_input_type), intent(in), dimension(:) :: input_comps
		integer, intent(out) :: Nlum, Ndust
		integer :: i
		Nlum = 0; Ndust = 0
		do i=1,size(input_comps, 1)
			select case(trim(input_comps(i)%comp_type_name))
			case("stellar"); Nlum = Nlum + 1
			case("dust"); Ndust = Ndust + 1
			case default
				 print*, "sellist tyypi pole olemas veel", trim(input_comps(i)%comp_type_name)
				 stop
			end select
		end do
	end subroutine loenda_komponentide_tyypide_arvud

end module all_comp_module