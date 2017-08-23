module all_comp_module
! 	use constants
	use comp_module
	
	type all_comp_type
		class(comp_type), allocatable, dimension(:) 			:: comp
		integer, allocatable, dimension(:) 					:: comp_im_ref !komponentide pildid... 2D massiiv kuna pilte voib olla mitut tyypi
		logical, allocatable, dimension(:)  				:: recalc_comp_images
		integer 											:: N_comp !kogu komponentide arv
		!accounting .. 											
		integer 											:: N_dust=0
		integer, dimension(:), allocatable 					:: dust_comp_id
		integer 											:: N_DM=0
		integer, dimension(:), allocatable 					:: DM_comp_id
		integer 											:: N_stellar=0
		integer, dimension(:), allocatable 					:: stellar_comp_id
		integer 											:: N_gas=0 
		integer, dimension(:), allocatable					:: gas_comp_id
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
		type(comp_input_type), intent(in), dimension(:), allocatable, target :: input_comps
		type(all_comp_type), intent(out) :: all_comp
		type(prof_par_list_type), pointer ::  par_list
		integer :: i

		all_comp%N_comp = size(input_comps, 1)
		allocate(all_comp%comp(1:all_comp%N_comp))
		do i=1,all_comp%N_comp
			all_comp%comp(i)%comp_prof_name = input_comps(i)%comp_prof_name
			all_comp%comp(i)%comp_name = input_comps(i)%comp_name
			all_comp%comp(i)%comp_type_name  = input_comps(i)%comp_type_name 
			all_comp%comp(i)%population_name = input_comps(i)%population_name
			all_comp%comp(i)%incl     = input_comps(i)%incl%val    
			all_comp%comp(i)%dist     = input_comps(i)%dist%val    
			all_comp%comp(i)%cnt_x    = input_comps(i)%cnt_x%val   
			all_comp%comp(i)%cnt_y    = input_comps(i)%cnt_y%val   
			all_comp%comp(i)%pos      = input_comps(i)%pos%val     
			all_comp%comp(i)%theta0   = input_comps(i)%theta0%val
			!
			!profiili allokeerimine ja parameetrite paika panemine... siia peab iga profiili korral tulema kontroll
			!
			if(trim(all_comp%comp(i)%comp_prof_name) == "Einasto") allocate(prof_Einasto_type::all_comp%comp(i)%prof_den)
			!
			!
			!
			par_list=>input_comps(i)%prof_pars
			do while(par_list%filled)
				call all_comp%comp(i)%prof_den%set_val(trim(par_list%par_name), par_list%par%val)
				if(associated(par_list%next)) then
					par_list => par_list%next
				else
					exit
				end if
			end do
! 			call all_comp%comp(i)%prof_den%init_profile()
! 			call init_comp(all_comp%comp(i))
		end do
	end subroutine convert_input_comp_to_all_comp
	
	subroutine asenda_viited(input_comps, all_comp)
		implicit none
		type(comp_input_type), intent(in), dimension(:), allocatable, target :: input_comps
		type(all_comp_type), intent(inout) :: all_comp
		type(prof_par_list_type), pointer ::  par_list
		integer :: i
		real(rk) :: tmp

		all_comp%N_comp = size(input_comps, 1)
		do i=1,all_comp%N_comp
			all_comp%comp(i)%incl     = all_comp%comp(input_comps(i)%incl%ref)%incl
			all_comp%comp(i)%dist      = all_comp%comp(input_comps(i)%dist%ref)%dist   
			all_comp%comp(i)%cnt_x     = all_comp%comp(input_comps(i)%cnt_x%ref)%cnt_x  
			all_comp%comp(i)%cnt_y     = all_comp%comp(input_comps(i)%cnt_y%ref)%cnt_y  
			all_comp%comp(i)%pos       = all_comp%comp(input_comps(i)%pos%ref)%pos    
			all_comp%comp(i)%theta0     = all_comp%comp(input_comps(i)%theta0%ref)%theta0
			par_list=>input_comps(i)%prof_pars
			do while(par_list%filled)
				call all_comp%comp( par_list%par%ref )%prof_den%get_val(trim(par_list%par_name), tmp) !votab v22rtuse ref jaoks
				call all_comp%comp(i)%prof_den%set_val(trim(par_list%par_name), tmp) !paneb v22rtuse teise kohta paika
				if(associated(par_list%next)) then
					par_list => par_list%next
				else
					exit
				end if
			end do
			call all_comp%comp(i)%prof_den%init_profile()
			call init_comp(all_comp%comp(i))
		end do
	end subroutine asenda_viited
	
	
end module all_comp_module