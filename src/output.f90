module output_module
	use likelihood_module


contains
	subroutine output_like_input(input_comps)
		implicit none
		integer :: mitmes_cube, i
		type(comp_input_type), dimension(:), allocatable, intent(in), target :: input_comps
		type(prof_par_list_type), pointer ::  par_list
		character(len=default_character_length) :: nimi
		integer :: iunit
		interface
			subroutine tryki_smth(nimi,par)
				import par_type_real
				import default_character_length
				implicit none
				character(len=default_character_length), intent(in) :: nimi
				type(par_type_real), intent(in) :: par
			end subroutine tryki_smth
		end interface
		procedure(tryki_smth), pointer :: tryki_output_par 
		
		tryki_output_par => tryki_prioriga
		if(kas_output_reana) tryki_output_par => tryki_reana
		
		iunit = 19
		open(file=output_fit_file, action="write", unit = iunit)
		do i=1,size(input_comps)
			if(.not.kas_output_reana) write(unit=iunit, fmt = "(A1,A,A1)")  "[",trim(input_comps(i)%comp_name),"]"
			if(.not.kas_output_reana) write(unit=iunit, fmt = "(A,A)") "prof = ", trim(input_comps(i)%comp_prof_name)
			if(.not.kas_output_reana) write(unit=iunit, fmt = "(A,A)") "type = ", trim(input_comps(i)%comp_type_name)
				nimi = "dist"
				call tryki_output_par(nimi, input_comps(i)%dist)
				nimi = "incl"
				call tryki_output_par(nimi, input_comps(i)%incl)
				nimi = "cnt_x"
				call tryki_output_par(nimi, input_comps(i)%cnt_x)
				nimi = "cnt_y"
				call tryki_output_par(nimi, input_comps(i)%cnt_y)
				nimi = "pos_wrt_phys_coord"
				call tryki_output_par(nimi, input_comps(i)%pos)
				nimi = "theta0"
				call tryki_output_par(nimi, input_comps(i)%theta0)
				par_list=>input_comps(i)%prof_pars
				do while(par_list%filled)
					call tryki_output_par(par_list%par_name, par_list%par)
					if(associated(par_list%next)) then
						par_list => par_list%next
					else; exit; end if
				end do
				if(.not.kas_output_reana) write(unit=iunit, fmt = "(A)") ""
		end do
		write(unit=iunit, fmt=*) ""
		close(iunit)
		nullify(par_list)
	contains
		subroutine tryki_prioriga(nimi,par)
			character(len=default_character_length), intent(in) :: nimi
			type(par_type_real), intent(in) :: par
			real(rk) :: kordaja
			select case(trim(nimi))
			case("incl"); kordaja = 180.0/pi
			case("cnt_x"); kordaja = 1.0/arcsec_to_rad
			case("cnt_y"); kordaja = 1.0/arcsec_to_rad
			case default; kordaja = 1.0
			end select
			write(unit=iunit, fmt = "(A22,A,E13.5,L,2E13.5,A1,A)") trim(nimi), " = ", par%val*kordaja, par%kas_fitib, par%min*kordaja, par%max*kordaja, " ", trim(input_comps(par%ref)%comp_name)
		end subroutine tryki_prioriga
		subroutine tryki_reana(nimi,par)
			character(len=default_character_length), intent(in) :: nimi
			type(par_type_real), intent(in) :: par
			real(rk) :: kordaja
			select case(trim(nimi))
			case("incl"); kordaja = 180.0/pi
			case("cnt_x"); kordaja = 1.0/arcsec_to_rad
			case("cnt_y"); kordaja = 1.0/arcsec_to_rad
			case default; kordaja = 1.0
			end select
			if(par%kas_fitib) write(unit=iunit, fmt="(E13.5,A)", advance="no") par%val*kordaja, " "
		end subroutine tryki_reana
	end subroutine output_like_input


	subroutine output_ML(input_comps, images)
		implicit none
		type(comp_input_type), dimension(:), allocatable, intent(inout), target :: input_comps
		type(image_type), dimension(:), allocatable, intent(in) :: images
		integer :: i,j, iunit
		iunit = 28
		open(file = output_ML_file, unit = iunit, action = "write")
		 
		write(unit=iunit, fmt="(A15)", advance = "no"), "Image"
		do j=1,size(input_comps)
			write(unit=iunit, fmt="(A15,A1)", advance = "no"), trim(input_comps(j)%comp_name), " "
		end do
		write(unit=iunit, fmt="(A15)", advance = "yes"), ""
		do i=1,size(images)
			write(unit=iunit, fmt="(A15)", advance = "no"), trim(images(i)%name)
			do j=1,size(input_comps)
				write(unit=iunit, fmt="(E15.8,A1)", advance = "no"), last_ML(i,j), " "
			end do
			write(unit=iunit, fmt="(A15)", advance = "yes"), ""
		end do
		close(unit=iunit)
	end subroutine output_ML
	
	subroutine output_images(input_comps, images)
		implicit none
		type(comp_input_type), dimension(:), allocatable, intent(inout), target :: input_comps
		type(image_type), dimension(:), allocatable, intent(in) :: images
		type(all_comp_type) :: all_comp
		real(rk) :: LL
		real(rk), dimension(:,:,:), allocatable :: v2ljundpildid
		real(rk), dimension(:,:), allocatable :: pilt
		integer :: i,j
		
		call convert_input_comp_to_all_comp(input_comps, all_comp)
		call asenda_viited(input_comps, all_comp) !all_comp muutujas asendamine

		select case(mis_fittimise_tyyp)
		case(2) !ehk componendid ja ML
			LL =  calc_log_likelihood(all_comp, images, v2ljundpildid)
			if(kas_koik_pildid_samast_vaatlusest) then
				allocate(pilt(1:size(v2ljundpildid, 2), 1:size(v2ljundpildid, 3)))
			else
				stop "Not implemented: v2ljundit ei saa teha kui koik pole sama suured pildid"
			end if
			do i=1,size(images)
				call write_matrix_to_fits(v2ljundpildid(i,:,:), images(i)%output_mdl_file)
				pilt = images(i)%obs - v2ljundpildid(i,:,:); where(.not.images(i)%mask) pilt = 0.0_rk !selleks, et maskeeriud t2hed ei s2raks hirmsasti
				call write_matrix_to_fits(pilt, images(i)%output_diff_file)
				pilt = abs(images(i)%obs - v2ljundpildid(i,:,:))/images(i)%sigma ; where(.not.images(i)%mask) pilt = 0.0_rk
				call write_matrix_to_fits( pilt, images(i)%output_rel_diff_file)
			end do
			if(.not.kas_vaikselt) print*, "Final LL was", LL
		case default
			print*, "no output available for this type", mis_fittimise_tyyp
		end select
		
	end subroutine output_images
end module output_module
	
	