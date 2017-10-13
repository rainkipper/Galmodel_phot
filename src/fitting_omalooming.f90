module fitting_omalooming_module
	use omalooming_module
	use likelihood_module
contains
	subroutine fittimine_omalooming(images, input_comps, all_comp)
		implicit none
		type(comp_input_type), dimension(:), allocatable, intent(inout), target :: input_comps
		type(all_comp_type), intent(out) :: all_comp !ehk v2ljundiks
		type(image_type), dimension(:), allocatable, intent(in) :: images
		real(rk), dimension(:), allocatable :: res, alumine, ylemine
		integer :: i, N
		type(prof_par_list_type), pointer ::  par_list
		integer mitmes_cube
		
		N = leia_vabade_parameetrite_arv(input_comps) !moodulist comp.f90
		allocate(res(1:N)); allocate(alumine(1:N)); allocate(ylemine(1:N))
		res = 0.0; alumine = 0.0; ylemine = 0.0;
		mitmes_cube = 0
		do i=1,size(input_comps)
			if(input_comps(i)%incl%kas_fitib)  then
				mitmes_cube=mitmes_cube+1
				alumine(mitmes_cube) = input_comps(i)%incl%min
				ylemine(mitmes_cube) = input_comps(i)%incl%max
			end if	
			if(input_comps(i)%cnt_x%kas_fitib) then
				mitmes_cube=mitmes_cube+1
				alumine(mitmes_cube) = input_comps(i)%cnt_x%min
				ylemine(mitmes_cube) = input_comps(i)%cnt_x%max	
			end if	
			if(input_comps(i)%cnt_y%kas_fitib)  then
				mitmes_cube=mitmes_cube+1
				alumine(mitmes_cube) = input_comps(i)%cnt_y%min
				ylemine(mitmes_cube) = input_comps(i)%cnt_y%max
			end if	
			if(input_comps(i)%pos%kas_fitib)  then
				mitmes_cube=mitmes_cube+1
				alumine(mitmes_cube) = input_comps(i)%pos%min
				ylemine(mitmes_cube) = input_comps(i)%pos%max
			end if	
			if(input_comps(i)%theta0%kas_fitib)  then
				mitmes_cube=mitmes_cube+1
				alumine(mitmes_cube) = input_comps(i)%theta0%min
				ylemine(mitmes_cube) = input_comps(i)%theta0%max
			end if	
			par_list=>input_comps(i)%prof_pars
			do while(par_list%filled)
				if(par_list%par%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					alumine(mitmes_cube) = par_list%par%min
					ylemine(mitmes_cube) = par_list%par%max
				end if	
				if(associated(par_list%next)) then
					par_list => par_list%next
				else
					exit
				end if
			end do
		end do
		nullify(par_list)
		
		
		res = fittimine(calc_LL, alumine, ylemine)
		print*, "omalooming lopetas"
		alumine(1) = calc_LL(res) !et saada input_comps paika
		
		
	contains
		function calc_LL(cube) result(res)
			implicit none
			real(rk) :: res
			real(rk), dimension(:), intent(in) :: Cube
			integer :: i
			type(prof_par_list_type), pointer ::  par_list
			integer mitmes_cube

			mitmes_cube = 0
			do i=1,size(input_comps)
				if(input_comps(i)%incl%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%incl%val = cube(mitmes_cube)
				end if	
				if(input_comps(i)%cnt_x%kas_fitib) then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%cnt_x%val = cube(mitmes_cube)	
				end if	
				if(input_comps(i)%cnt_y%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%cnt_y%val = cube(mitmes_cube)
				end if	
				if(input_comps(i)%pos%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%pos%val = cube(mitmes_cube)
				end if	
				if(input_comps(i)%theta0%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%theta0%val = cube(mitmes_cube)
				end if	
				par_list=>input_comps(i)%prof_pars
				do while(par_list%filled)
					if(par_list%par%kas_fitib)  then
						mitmes_cube=mitmes_cube+1
						par_list%par%val = cube(mitmes_cube)
					end if	
					if(associated(par_list%next)) then
						par_list => par_list%next
					else
						exit
					end if
				end do
			end do
			nullify(par_list)

			call convert_input_comp_to_all_comp(input_comps, all_comp)
			call asenda_viited(input_comps, all_comp) !all_comp muutujas asendamine
			
! 			call prindi_koik_fititavad_parameetrid(input_comps)
			
			res = calc_log_likelihood(all_comp, images) 
			!kui lisamassidele juurde asju arvutatud, siis paneb uued massid vastavalt eelmistele... input_comps juurde
			nullify(par_list)
		end function calc_LL
	end subroutine fittimine_omalooming
end module fitting_omalooming_module