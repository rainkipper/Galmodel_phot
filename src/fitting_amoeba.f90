module fitting_amoeba_module
	use mynr_amoeba_module
	use likelihood_module
	private
	public fittimine_amoeba
contains
	subroutine fittimine_amoeba(images, input_comps, all_comp, input_points)
		implicit none
		!
		! jagatud sisend-v2ljund teiste fittijatega
		!
		type(comp_input_type), dimension(:), allocatable, intent(inout), target :: input_comps
		type(all_comp_type), intent(inout) :: all_comp !ehk v2ljundiks
		type(image_type), dimension(:), allocatable, intent(in) :: images
		real(rk), intent(in), optional, dimension(:,:), allocatable :: input_points !valikuline sisendpunktide valik... juhul kui tahta multinestist edasi teha	
		real(rk), dimension(:,:), allocatable :: points
		real(rk), dimension(:), allocatable :: likelihoods 
		integer :: N_vabu_par, j, iter
		
		
		
		N_vabu_par = leia_vabade_parameetrite_arv(input_comps)
		allocate(points(N_vabu_par+1, N_vabu_par)); allocate(likelihoods(N_vabu_par))
		if(present(input_points)) then
			if(size(input_points,1).ne.size(input_points,1) .or. size(input_points,2).ne.size(input_points,2)) then
				stop "Amoeba juures sisendiks antud punktide suurus ei klapi vajadusega"
			end if
			points = input_points
		else
			points = suvalised_priori_punktid(input_comps)
		end if
		
		do j=1,size(points,1)
			likelihoods(j) = calc_LL(points(j,:))
		end do
		print*, "to amoeba"
		call my_amoeba(points, likelihoods, amoeba_fractional_tolerance, calc_LL, iter)
		
	contains
		function suvalised_priori_punktid(input_comps) result(points)
			implicit none
			type(comp_input_type), dimension(:), allocatable, intent(in), target :: input_comps
			real(rk), dimension(:,:), allocatable :: points
			real(rk), dimension(:), allocatable :: single_par_random
			integer :: i, mitmes_cube
			type(prof_par_list_type), pointer ::  par_list
			allocate(single_par_random(1:(N_vabu_par+1)))
			if(.not.allocated(points)) allocate(points(1:(N_vabu_par+1), 1:N_vabu_par))
			mitmes_cube = 0
			do i=1,size(input_comps)
				if(input_comps(i)%incl%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					call random_number(single_par_random)
					points(:,mitmes_cube) = get_random_to_prior(single_par_random, input_comps(i)%incl)
				end if	
				if(input_comps(i)%cnt_x%kas_fitib) then
					mitmes_cube=mitmes_cube+1
					call random_number(single_par_random)
					points(:,mitmes_cube) = get_random_to_prior(single_par_random, input_comps(i)%cnt_x)
				end if	
				if(input_comps(i)%cnt_y%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					call random_number(single_par_random)
					points(:,mitmes_cube) = get_random_to_prior(single_par_random, input_comps(i)%cnt_y)
				end if	
				if(input_comps(i)%pos%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					call random_number(single_par_random)
					points(:,mitmes_cube) = get_random_to_prior(single_par_random, input_comps(i)%pos)
				end if	
				if(input_comps(i)%theta0%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					call random_number(single_par_random)
					points(:,mitmes_cube) = get_random_to_prior(single_par_random, input_comps(i)%theta0)
				end if	
				par_list=>input_comps(i)%prof_pars
				do while(par_list%filled)
					if(par_list%par%kas_fitib)  then
						mitmes_cube=mitmes_cube+1
						call random_number(single_par_random)
						points(:,mitmes_cube) = get_random_to_prior(single_par_random, par_list%par)
					end if	
						if(associated(par_list%next)) then
							par_list => par_list%next
						else
							exit
					end if
				end do
			end do
			nullify(par_list)
		end function suvalised_priori_punktid
		elemental function get_random_to_prior(rnd, prior) result(res)
			implicit none
			real(rk), intent(in) :: rnd
			type(par_type_real), intent(in) :: prior
			real(rk) :: res
			real(rk), parameter :: kui_keskelt = 0.2
! 			res = rnd*(prior%max - prior%min) + prior%min !yhtlaselt yle terve priori
			res = (rnd-0.5)*(prior%max - prior%min)*kui_keskelt + 0.5*(prior%min+prior%max) !priori keskemalt
		end function get_random_to_prior
		function calc_LL(cube) result(res)
			implicit none
			real(rk) :: res
			real(rk), dimension(:), intent(in) :: Cube
			integer :: i
			type(prof_par_list_type), pointer ::  par_list
			integer mitmes_cube
			real(rk), dimension(:), allocatable :: lisakaalud_massile

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
			
			res = -1.0* calc_log_likelihood(all_comp, images) !-1* kuna leiab miinimumi
			!kui lisamassidele juurde asju arvutatud, siis paneb uued massid vastavalt eelmistele... input_comps juurde
			nullify(par_list)
		end function calc_LL
	end subroutine fittimine_amoeba
	subroutine prindi_koik_fititavad_parameetrid(input_comps)
		implicit none
		integer :: mitmes_cube, i
		type(comp_input_type), dimension(:), allocatable, intent(in), target :: input_comps
		type(prof_par_list_type), pointer ::  par_list
		character(len=default_character_length) :: nimi
			
		mitmes_cube = 0
		do i=1,size(input_comps)
			if(input_comps(i)%incl%kas_fitib)  then
				mitmes_cube=mitmes_cube+1
				nimi = "incl"
				call tryki_yks(nimi, input_comps(i)%incl)
			end if	
			if(input_comps(i)%cnt_x%kas_fitib) then
				mitmes_cube=mitmes_cube+1
				nimi = "cntx"
				call tryki_yks(nimi, input_comps(i)%cnt_x)
! 				input_comps(i)%cnt_x%val = cube(mitmes_cube)
			end if	
			if(input_comps(i)%cnt_y%kas_fitib)  then
				mitmes_cube=mitmes_cube+1
				nimi = "cnty"
				call tryki_yks(nimi, input_comps(i)%cnt_y)
			end if	
			if(input_comps(i)%pos%kas_fitib)  then
				mitmes_cube=mitmes_cube+1
				nimi = "pos"
				call tryki_yks(nimi, input_comps(i)%pos)
			end if	
			if(input_comps(i)%theta0%kas_fitib)  then
				mitmes_cube=mitmes_cube+1
				nimi = "theta0"
				call tryki_yks(nimi, input_comps(i)%theta0)
			end if	
			par_list=>input_comps(i)%prof_pars
			do while(par_list%filled)
				if(par_list%par%kas_fitib)  then
					mitmes_cube=mitmes_cube+1					
					call tryki_yks(par_list%par_name, par_list%par)
				end if	
				if(associated(par_list%next)) then
					par_list => par_list%next
				else
					exit
				end if
			end do
		end do
		nullify(par_list)
	contains
		subroutine tryki_yks(nimi,par)
			character(len=default_character_length), intent(in) :: nimi
			type(par_type_real), intent(in) :: par
			print "(A,3(A,F10.5))", trim(nimi), ": ", par%min, " ----- ", par%val, " ------ ", par%max
		end subroutine tryki_yks
	end subroutine prindi_koik_fititavad_parameetrid
end module fitting_amoeba_module