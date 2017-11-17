module read_input_module
	!
	!TODO vaja veel lisada accounting osa all_comp sisse lugemisel
	!TODO iga faili sisselugemiel vaja teha kontroll, kas koik olulised asjad on sisse loetud.
	!TODO praegu on siin testfiltrid... vaja oleks need asjendada p2ris filtrite sisselugemisega
	!
	use all_comp_module
	use filters_module
	use images_module
	use file_operations_module
	use populations_module
	use filters_module
	use dust_module
	use filters_module

	character(len=default_character_length), dimension(:), allocatable :: viitamiseks_comp_nimed !vastavalt selline j2rjekord, mis input_comp failis on... globaalne muutuja yldkasutuseks eri asjade sisselugemisel
! 	type(prof_par_list_type), dimension(1:max_prof_par_size), target :: all_prof_pars
! 	integer :: mitu_all_prof_pars_kasutatud = 0
contains
	subroutine read_components(comp_file, input_comps) 
		implicit none
		type(comp_input_type), intent(out), dimension(:), allocatable, target :: input_comps
		character(len=default_character_length), intent(in) :: comp_file
		integer :: N_comp
		integer :: id1, id2, mitmes_comp = 0
		integer :: ios
		integer, parameter :: iunit = 11
		character(len=default_character_length) :: rida, subrida
		type(prof_par_list_type), pointer ::  par_list
		
		!
		N_comp = read_nr_of_comp(comp_file);
		allocate(input_comps(1:N_comp))
		call read_list_of_references(); !print*, viitamiseks_comp_nimed
		

		open(file = comp_file, unit = iunit, action = "read")
		ios = 1
		do while(ios .ge. 0)
			rida = repeat(" ", default_character_length) !m2lu nullimine ymberkasutamiseks
			read(fmt = "(A)", iostat=ios, unit=iunit) rida
			if( len(trim(rida))<2 .or. rida(1:1)=="#") cycle !mottetute ridade v2ltimine
			if(rida(1:1) == "[") then
				mitmes_comp = mitmes_comp + 1
				input_comps(mitmes_comp)%comp_name = rida(2:(index(rida, "]")-1))
				continue
			end if
			id1 = 1; id2 = index(rida, "=")
			if(id2==0) cycle
! 			print*, "suundub ---> ", trim(rida(1:(id2-1)))
			subrida = repeat(" ", default_character_length)
			subrida = trim(ADJUSTL(rida((id2+1):len(rida))))
			select case(trim(rida(1:(id2-1))))
			case("prof")
				input_comps(mitmes_comp)%comp_prof_name = subrida
			case("type")
				input_comps(mitmes_comp)%comp_type_name = subrida
			case("population")
				input_comps(mitmes_comp)%population_name = subrida
			case("dist")
				input_comps(mitmes_comp)%dist = read_real_par(subrida)
			case("theta0")
				input_comps(mitmes_comp)%theta0 = read_real_par(subrida)
				input_comps(mitmes_comp)%theta0%val = input_comps(mitmes_comp)%theta0%val * pi/180.0
				input_comps(mitmes_comp)%theta0%min = input_comps(mitmes_comp)%theta0%min * pi/180.0
				input_comps(mitmes_comp)%theta0%max = input_comps(mitmes_comp)%theta0%max * pi/180.0
			case("incl")
				input_comps(mitmes_comp)%incl = read_real_par(subrida)
				!vaja konvertida ka radiaanidesse
				input_comps(mitmes_comp)%incl%val = input_comps(mitmes_comp)%incl%val * pi/180.0
				input_comps(mitmes_comp)%incl%min = input_comps(mitmes_comp)%incl%min * pi/180.0
				input_comps(mitmes_comp)%incl%max = input_comps(mitmes_comp)%incl%max * pi/180.0
			case("pos_wrt_phys_coord")
				input_comps(mitmes_comp)%pos = read_real_par(subrida)
				input_comps(mitmes_comp)%pos%val = input_comps(mitmes_comp)%pos%val * pi/180.0
				input_comps(mitmes_comp)%pos%min = input_comps(mitmes_comp)%pos%min * pi/180.0
				input_comps(mitmes_comp)%pos%max = input_comps(mitmes_comp)%pos%max * pi/180.0
			case("cnt_x")
				input_comps(mitmes_comp)%cnt_x = read_real_par(subrida)				 !see on olemuselt nurk ning radiaanides
				input_comps(mitmes_comp)%cnt_x%val = arcsec_to_rad * input_comps(mitmes_comp)%cnt_x%val
				input_comps(mitmes_comp)%cnt_x%min = arcsec_to_rad * input_comps(mitmes_comp)%cnt_x%min
				input_comps(mitmes_comp)%cnt_x%max = arcsec_to_rad * input_comps(mitmes_comp)%cnt_x%max

			case("cnt_y")
				input_comps(mitmes_comp)%cnt_y = read_real_par(subrida)
				input_comps(mitmes_comp)%cnt_y%val = arcsec_to_rad * input_comps(mitmes_comp)%cnt_y%val
				input_comps(mitmes_comp)%cnt_y%min = arcsec_to_rad * input_comps(mitmes_comp)%cnt_y%min
				input_comps(mitmes_comp)%cnt_y%max = arcsec_to_rad * input_comps(mitmes_comp)%cnt_y%max
			case default
				par_list => input_comps(mitmes_comp)%prof_pars
				do while(associated(par_list%next))
					par_list => par_list%next
				end do
				if(par_list%filled) then
					allocate(par_list%next)
					par_list=>par_list%next
				end if
				par_list%par_name = trim(rida(1:(id2-1)))
				par_list%par = read_real_par(subrida)
! print "(A,A,A,A,I)", trim(input_comps(mitmes_comp)%comp_name), ":", trim(par_list%par_name), "->",par_list%par%ref
				par_list%filled = .true.
				nullify(par_list%next)
			end select
		end do
		close(unit = iunit)

	contains
		subroutine read_list_of_references() 
			!eesm2rk on muuta sonalised komponentide viited (nagu "bulge") numbrilisteks komponentide numbriteks
			implicit none
			integer :: id1, id2
			integer :: ios
			integer, parameter :: iunit = 10
			integer :: mitmes_comp
			mitmes_comp = 0
			allocate(viitamiseks_comp_nimed(1:N_comp)) !algne pealmisest subroutine-st
			viitamiseks_comp_nimed(:) = repeat(" ", default_character_length)
			open(file = comp_file, unit = iunit, action = "read")
			do while(ios .ge. 0)
				rida = repeat(" ", default_character_length) !m2lu nullimine ymberkasutamiseks
				read(fmt = "(A)", iostat=ios, unit=iunit) rida
				if(rida(1:1) == "[") then
					mitmes_comp = mitmes_comp + 1
					viitamiseks_comp_nimed(mitmes_comp) = trim(rida(2:(index(rida, "]")-1)))
				end if
			end do
			close(unit = iunit)
		end subroutine read_list_of_references
		function get_ref_from_name(vidin) result(res)
			implicit none
			character(len=default_character_length), intent(in) :: vidin
			integer :: res, i
			res = -1
			do i=1,N_comp
				if(trim(vidin) == trim(viitamiseks_comp_nimed(i))) res = i
			end do
			if(res == -1) stop "Niisugust comp nime pole olemas faili sisselugemisel (veel)"
		end function get_ref_from_name
		
		function read_real_par(rida) result(res)
			implicit none
			character(len=default_character_length), intent(in) :: rida
			character(len=default_character_length) :: pisirida
			type(par_type_real) :: res
			pisirida = repeat(" ", default_character_length) !m2lu nullimiseks
			read(rida, fmt=*) res%val, res%kas_fitib, res%min, res%max, pisirida
			res%ref = get_ref_from_name(pisirida)
		end function read_real_par
	end subroutine read_components
	subroutine read_images_and_filters(images_file, images)
		implicit none
		type(image_type), intent(out), dimension(:), allocatable :: images
		character(len=default_character_length), intent(in) :: images_file
! 		type(filter_type), dimension(:), allocatable, intent(out), target :: filters !vaja, et linkida kohe kylge
		integer :: N_im
		integer :: i
		integer :: id1, id2 !ajutised abimuutujad	
		integer :: mitmes_im
		integer :: iunit, ios !molemad failindusega seotud muutujad
		character(len=default_character_length) :: rida, subrida, tyhi
		real(rk), dimension(:,:), allocatable :: tmp_pilt !ehk teha loogilist pilti
		

		N_im = read_nr_of_comp(images_file)
		allocate(images(1:N_im))
		open(file = images_file, unit = iunit, action = "read")
		ios = 1
		do while(ios .ge. 0)
			rida = tyhi !m2lu nullimine ymberkasutamiseks... selleks, et trim k2suga asjade allaj22nud asju ei kasutaks
			read(fmt = "(a)", iostat=ios, unit=iunit) rida
			rida = adjustl(rida) !eemdalab koik esimesed tyhikud ja paneb loppu... 
			if( len(trim(rida))<2 .or. rida(1:1)=="#") cycle !mottetute ridade v2ltimine
			if(rida(1:1) == "[") then
				mitmes_im = mitmes_im + 1
				images(mitmes_im)%name = tyhi
				images(mitmes_im)%name = trim(adjustl(rida(2:(index(rida, "]")-1))))
				continue
			end if
			id1 = 1; id2 = index(rida, "=")
			if(id2==0) cycle
! 			print*, "suundub ---> ", trim(rida(1:(id2-1)))
			subrida = tyhi
			subrida = trim(ADJUSTL(rida((id2+1):len(rida))))
			select case(trim(adjustl(rida(1:(id2-1)))))
			case("filter")
				nullify(images(mitmes_im)%filter)
				do i=1, size(filters)
					if(trim(filters(i)%name)==trim(subrida)) then
						images(mitmes_im)%filter => filters(i)
						exit
					end if
				end do
				if(.not.associated(images(mitmes_im)%filter)) then
					print*, "Niisugust pilti ei saanud mingi filtriga vastavusse viia: ", trim(subrida)
					stop  !kontroll, kas ikkagi leidis mingi filtri vaste
				end if
			case("scale_x"); read(subrida, fmt=*) images(mitmes_im)%scale_x
			case("sky_noise"); read(subrida, fmt=*) images(mitmes_im)%sky_noise
			case("scale_y"); read(subrida, fmt=*) images(mitmes_im)%scale_y
			case("angle_wrt_phys_coord")
				read(subrida, fmt=*) images(mitmes_im)%pos
				images(mitmes_im)%pos = images(mitmes_im)%pos * pi / 180.0
				images(mitmes_im)%cos_pos = cos(images(mitmes_im)%pos) !vaja ainult yhe korra arvutada
				images(mitmes_im)%sin_pos = sin(images(mitmes_im)%pos) !vaja ainult yhe korra arvutada
				images(mitmes_im)%tan_pos = tan(images(mitmes_im)%pos) !vaja ainult yhe korra arvutada
			case("image_coord_cnt_x"); read(subrida, fmt=*) images(mitmes_im)%x
			case("image_coord_cnt_y"); read(subrida, fmt=*) images(mitmes_im)%y
			case("obs_file"); 
				read(subrida, fmt="(a)") images(mitmes_im)%obs_file
				call read_fits_to_matrix(trim(images(mitmes_im)%obs_file), images(mitmes_im)%obs)
			case("sigma_file"); 
				read(subrida, fmt="(a)") images(mitmes_im)%sigma_file
				call read_fits_to_matrix(trim(images(mitmes_im)%sigma_file), images(mitmes_im)%sigma)	
			case("mask_file"); 
				read(subrida, fmt="(a)") images(mitmes_im)%mask_file
				call read_fits_to_matrix(trim(images(mitmes_im)%mask_file), tmp_pilt)
				allocate(images(mitmes_im)%mask(1:size(tmp_pilt,1), 1:size(tmp_pilt,2)))
				images(mitmes_im)%mask = tmp_pilt>0.5_rk !ehk yle 0.5 kasutatakse
				deallocate(tmp_pilt)
			case("psf_file"); 
				read(subrida, fmt="(a)") images(mitmes_im)%psf_file
				call read_fits_to_matrix(trim(images(mitmes_im)%psf_file), images(mitmes_im)%psf)	 
			case("mdl_file"); 
				read(subrida, fmt="(a)") images(mitmes_im)%output_mdl_file
			case("diff_file"); 
				read(subrida, fmt="(a)") images(mitmes_im)%output_diff_file
			case("rel_diff_file"); 
				read(subrida, fmt="(a)") images(mitmes_im)%output_rel_diff_file
			case default
				print*, "Niisugust asja ei saa piltide juures lugeda: ",trim(adjustl(rida(1:(id2-1))))
				stop
			end select
		end do
		close(unit = iunit)		
	end subroutine read_images_and_filters
	function read_nr_of_comp(mis_failist) result(res)
		implicit none
		integer, parameter :: iunit = 12
		integer :: ios
		integer :: res
		character(len=default_character_length) :: rida
		character(len=default_character_length), intent(in) :: mis_failist
		ios  = 0
		res = 0
		open(file = trim(mis_failist), unit = iunit, action = "read")
		do while(ios .ge. 0)
			rida = repeat(" ", default_character_length) !m2lu nullimine ymberkasutamiseks
			read(fmt = "(A)", iostat=ios, unit=iunit) rida
			if(rida(1:1)=="[") res = res + 1
		end do
		close(unit = iunit)
	end function read_nr_of_comp

end module read_input_module

