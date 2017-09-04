program gm
	use fill_model_image_module
	use	read_input_module
	use fitting_multinest_module
	!  	use fitsio

	type(all_comp_type) :: all_comp, all_comp_input
	type(comp_input_type), dimension(:), allocatable :: input_comp
	type(model_image_real_type) :: mdl1, mdl2, mdl
	type(image_type), dimension(:), allocatable :: images
	type(filter_type), dimension(:), allocatable :: filters
	integer, parameter :: N_comp = 1
	integer :: test_image_number
	character(len=*), parameter :: fname="comp_im.txt"
	character(len=default_character_length), parameter :: input_comp_file = "Input/input_comp.txt"
	character(len=default_character_length), parameter :: input_image_file = "Input/input_images.txt"
	logical :: kas_comp_im, kas_los
	integer :: i, N
	real(rk) :: t1,t2 !aja mootmine
	real(rk) :: ll
	logical, parameter :: kas_fitib_vs_lihtsalt_ouput = .false.
	
	call cpu_time(t1)
	!
	! ================== mingi testgalaktika tekitamine =================
	!
	call read_components(input_comp_file, input_comp)
	call read_images_and_filters(input_image_file, images, filters)
	
	if(kas_fitib_vs_lihtsalt_ouput) then
		call jooksuta_fittimine(images, input_comp, all_comp)
	else
		call convert_input_comp_to_all_comp(input_comp, all_comp)
		call asenda_viited(input_comp, all_comp)
		call create_model_image_from_obs(mdl, images(1))
		call write_matrix_to_fits(mdl%mx, images(1)%output_mdl_file)
	end if
	

! 	ll = calc_log_likelihood(all_comp, images)
	
	call  cpu_time(t2)
	print*, "All done:D", t2-t1

contains
    
end program gm