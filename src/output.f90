module output_module
	use likelihood_module
contains
	subroutine output_like_input(input_comps)
		implicit none
		integer :: mitmes_cube, i
		type(comp_input_type), dimension(:), allocatable, intent(in), target :: input_comps
		type(prof_par_list_type), pointer ::  par_list
		character(len=default_character_length) :: nimi
			
		do i=1,size(input_comps)
			print "(A1,A,A1)", "[",trim(input_comps(i)%comp_name),"]"
			print "(A11,A)", "prof = ", trim(input_comps(i)%comp_prof_name)
			print "(A11,A)", "type = ", trim(input_comps(i)%comp_type_name)
				nimi = "dist"
				call tryki_yks(nimi, input_comps(i)%dist)
				nimi = "incl"
				call tryki_yks(nimi, input_comps(i)%incl)
				nimi = "cnt_x"
				call tryki_yks(nimi, input_comps(i)%cnt_x)
				nimi = "cnt_y"
				call tryki_yks(nimi, input_comps(i)%cnt_y)
				nimi = "pos_wrt_phys_coord"
				call tryki_yks(nimi, input_comps(i)%pos)
! 				nimi = "theta0"
! 				call tryki_yks(nimi, input_comps(i)%theta0)
				par_list=>input_comps(i)%prof_pars
				do while(par_list%filled)
					call tryki_yks(par_list%par_name, par_list%par)
					if(associated(par_list%next)) then
						par_list => par_list%next
					else; exit; end if
				end do
				print*, ""
		end do
		nullify(par_list)
	contains
		subroutine tryki_yks(nimi,par)
			character(len=default_character_length), intent(in) :: nimi
			type(par_type_real), intent(in) :: par
			print "(A8,A,E13.5,L,2E13.5,A1,A)", trim(nimi), " = ", par%val, par%kas_fitib, par%min, par%max, " ", trim(input_comps(par%ref)%comp_name)
		end subroutine tryki_yks
	end subroutine output_like_input
	subroutine output_ML(input_comps, images)
		implicit none
		type(comp_input_type), dimension(:), allocatable, intent(inout), target :: input_comps
		type(image_type), dimension(:), allocatable, intent(in) :: images
		integer :: i,j
		
		write(unit=*, fmt="(A15)", advance = "no"), "Image"
		do j=1,size(input_comps)
			write(unit=*, fmt="(A15,A1)", advance = "no"), trim(input_comps(j)%comp_name), " "
		end do
		write(unit=*, fmt="(A15)", advance = "yes"), ""
		do i=1,size(images)
			write(unit=*, fmt="(A15)", advance = "no"), trim(images(i)%name)
			do j=1,size(input_comps)
				write(unit=*, fmt="(E15.8,A1)", advance = "no"), last_ML(i,j), " "
			end do
			write(unit=*, fmt="(A15)", advance = "yes"), ""
		end do
	end subroutine output_ML
end module output_module
	
	