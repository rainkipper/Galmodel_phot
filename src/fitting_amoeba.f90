module fitting_amoeba_module
	use mynr_amoeba_module
	use likelihood_module
	private
	public fittimine_amoeba
	real(rk), private, parameter :: ftol = 0.01
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
		call my_amoeba(points, likelihoods, ftol, calc_LL, iter)
		
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
					points(:,mitmes_cube) = single_par_random * (input_comps(i)%incl%max - input_comps(i)%incl%min) + input_comps(i)%incl%min
				end if	
				if(input_comps(i)%cnt_x%kas_fitib) then
					mitmes_cube=mitmes_cube+1
					call random_number(single_par_random)
					points(:,mitmes_cube) = single_par_random * (input_comps(i)%cnt_x%max - input_comps(i)%cnt_x%min) + input_comps(i)%cnt_x%min	
				end if	
				if(input_comps(i)%cnt_y%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					call random_number(single_par_random)
					points(:,mitmes_cube) = single_par_random * (input_comps(i)%cnt_y%max - input_comps(i)%cnt_y%min) + input_comps(i)%cnt_y%min
				end if	
				if(input_comps(i)%pos%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					call random_number(single_par_random)
					points(:,mitmes_cube) = single_par_random * (input_comps(i)%pos%max - input_comps(i)%pos%min) + input_comps(i)%pos%min
				end if	
				if(input_comps(i)%theta0%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					call random_number(single_par_random)
					points(:,mitmes_cube) = single_par_random * (input_comps(i)%theta0%max - input_comps(i)%theta0%min) + input_comps(i)%theta0%min
				end if	
				par_list=>input_comps(i)%prof_pars
				do while(par_list%filled)
					if(par_list%par%kas_fitib)  then
						mitmes_cube=mitmes_cube+1
						call random_number(single_par_random)
						points(:,mitmes_cube) = single_par_random * (par_list%par%max - par_list%par%min) + par_list%par%min
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
			res = -1.0* calc_log_likelihood(all_comp, images, lisakaalud_massile) !-1* kuna leiab miinimumi
			!kui lisamassidele juurde asju arvutatud, siis paneb uued massid vastavalt eelmistele... input_comps juurde
			do i=1,all_comp%N_comp
				par_list => input_comps(i)%prof_pars
				do while(par_list%filled)
					if(trim(par_list%par_name)=="M") then
						par_list%par%val = par_list%par%val * lisakaalud_massile(i)
						exit
					end if
					par_list => par_list%next
				end do
			end do
			nullify(par_list)
		end function calc_LL
	end subroutine fittimine_amoeba
end module fitting_amoeba_module