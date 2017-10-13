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
	
	subroutine output_images(input_comps, images)
		implicit none
		type(comp_input_type), dimension(:), allocatable, intent(inout), target :: input_comps
		type(image_type), dimension(:), allocatable, intent(in) :: images
		type(all_comp_type) :: all_comp
		real(rk) :: LL
		real(rk), dimension(:,:,:), allocatable :: v2ljundpildid
		integer :: i,j

		call convert_input_comp_to_all_comp(input_comps, all_comp)
		call asenda_viited(input_comps, all_comp) !all_comp muutujas asendamine
		
		select case(mis_fittimise_tyyp)
		case(2) !ehk componendid ja ML
			LL =  calc_log_likelihood(all_comp, images, v2ljundpildid)
			do i=1,size(images)
				call write_matrix_to_fits(v2ljundpildid[i,:,:], images(i)%output_mdl_file)
				call write_matrix_to_fits(images(i)%obs - v2ljundpildid[i,:,:], images(i)%output_diff_file)
				call write_matrix_to_fits( abs(images(i)%obs - v2ljundpildid[i,:,:])/images(i)%sigma , images(i)%output_rel_diff_file)
			end do
			print*, "Final LL was", LL
		case default
			print*, "no output available for this type"
		end select
		
	end subroutine output_images
end module output_module
	
	