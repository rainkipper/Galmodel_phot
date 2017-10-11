module all_comp_module
! 	use constants
	use comp_module
	
	type all_comp_type
		class(comp_type), allocatable, dimension(:) 			:: comp
		integer, allocatable, dimension(:) 					:: comp_im_ref !komponentide pildid... 2D massiiv kuna pilte voib olla mitut tyypi
		logical, allocatable, dimension(:)  				:: recalc_adaptive_ims
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
		type(all_comp_type), intent(inout) :: all_comp
		type(prof_par_list_type), pointer ::  par_list
		integer :: i

		all_comp%N_comp = size(input_comps, 1)
		
		if(not(allocated(all_comp%comp))) then
			allocate(all_comp%comp(1:all_comp%N_comp))
		end if
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
			if(trim(all_comp%comp(i)%comp_prof_name) == "Einasto" .and. not(allocated(all_comp%comp(i)%prof_den))) allocate(prof_Einasto_type::all_comp%comp(i)%prof_den)
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
			nullify(par_list)
! 			call all_comp%comp(i)%prof_den%init_profile()
! 			call init_comp(all_comp%comp(i))
		end do
	end subroutine convert_input_comp_to_all_comp
	
	subroutine asenda_viited(input_comps, all_comp, recalc_comp)
		implicit none
		type(comp_input_type), intent(in), dimension(:), allocatable, target :: input_comps
		type(all_comp_type), intent(inout) :: all_comp
		logical, intent(inout), dimension(:), allocatable, optional :: recalc_comp !siit lisatakse viidete kaudu ymberarvutuse vajadus
		type(prof_par_list_type), pointer ::  par_list
		integer :: i!,j
		integer :: viide
		real(rk) :: tmp!, tmp2
		logical :: kas_recalc !ehk kas on vaja t2ita ymber arvutuste asju
		kas_recalc = present(recalc_comp)
		all_comp%N_comp = size(input_comps, 1)
		do i=1,all_comp%N_comp
			
			viide = input_comps(i)%incl%ref
			all_comp%comp(i)%incl     = all_comp%comp(viide)%incl
! 			if(kas_recalc .and. input_comps(viide)%incl%kas_fitib .and.  recalc_comp(viide) ) recalc_comp(i) = .true.

			viide = input_comps(i)%incl%ref
			all_comp%comp(i)%dist      = all_comp%comp(viide)%dist
! 			if(kas_recalc .and. input_comps(viide)%dist%kas_fitib .and.  recalc_comp(viide) ) recalc_comp(i) = .true.
			   
			viide = input_comps(i)%cnt_x%ref
			all_comp%comp(i)%cnt_x     = all_comp%comp(viide)%cnt_x  
! 			if(kas_recalc .and. input_comps(viide)%cnt_x%kas_fitib .and.  recalc_comp(viide) ) recalc_comp(i) = .true.

			viide = input_comps(i)%cnt_y%ref
			all_comp%comp(i)%cnt_y     = all_comp%comp(viide)%cnt_y  
! 			if(kas_recalc .and. input_comps(viide)%cnt_y%kas_fitib .and.  recalc_comp(viide) ) recalc_comp(i) = .true.

			viide = input_comps(i)%pos%ref
			all_comp%comp(i)%pos       = all_comp%comp(viide)%pos    
! 			if(kas_recalc .and. input_comps(viide)%pos%kas_fitib .and.  recalc_comp(viide) ) recalc_comp(i) = .true.
			
			viide = input_comps(i)%theta0%ref
			all_comp%comp(i)%theta0     = all_comp%comp(viide)%theta0
! 			if(kas_recalc .and. input_comps(viide)%theta0%kas_fitib .and.  recalc_comp(viide) ) recalc_comp(i) = .true.

			par_list=>input_comps(i)%prof_pars
			do while(par_list%filled)
				call all_comp%comp( par_list%par%ref )%prof_den%get_val(trim(par_list%par_name), tmp) !votab v22rtuse ref jaoks
				call all_comp%comp(i)%prof_den%set_val(trim(par_list%par_name), tmp) !paneb v22rtuse teise kohta paika
				!ymber arvutuse kontrolli asjad
! 				viide = par_list%par%ref
! 				par_list2 => input_comps(viide)%prof_pars
! 				if(kas_recalc .and. recalc_comp(viide) ) then
! 					!vastava parameetri v2lja otsimine
! 					!seda konstruktsiooni vaja eristamaks  mittefititavat parameetrit fititavast
! 					do j=1,1000
! 						if( trim(par_list%par_name) == trim(par_list2%par_name) ) then
! 							exit
! 						else
! 							par_list2 => par_list2%next
! 						end if
! 					end do
! ! 					if(par_list2%par%kas_fitib) recalc_comp(i) = .true.
!
! 				end if
				!
				if(associated(par_list%next)) then
					par_list => par_list%next
				else
					exit
				end if
			end do
			call all_comp%comp(i)%prof_den%init_profile()
			call init_comp(all_comp%comp(i))
		end do
		nullify(par_list)
	end subroutine asenda_viited
	

	
end module all_comp_module