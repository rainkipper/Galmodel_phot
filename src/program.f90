program gm
	use fitting_module
	use read_input_module
	
	type(all_comp_type) :: all_comp
	type(comp_input_type), dimension(:), allocatable :: input_comp
	type(image_type), dimension(:), allocatable :: images
	integer, parameter :: N_comp = 1
	character(len=default_character_length), parameter :: input_comp_file = "input_comp.txt"
	character(len=default_character_length), parameter :: input_param_file = "input_parameters.txt"
	character(len=default_character_length), parameter :: input_image_file = "input_images.txt"
	real(rk) :: t1,t2 !aja mootmine


	call random_seed()
	call cpu_time(t1)
	call read_constants_and_parameters(input_param_file)
	call init_lum_dependencies(sol_file_name, populations_file_name) !siin s2titakse filtrid ja tolmu ja populatsioonide spektraaljaotused paika
	print*, "comp lugemisse"
	call read_components(input_comp_file, input_comp)
	print*, "pilte lugema"
	call read_images_and_filters(input_image_file, images)
	print*, "to init images"
	call init_images(images) !peamiselt psf crop		
	call fit_galaxy(images, input_comp, all_comp)
	call cpu_time(t2)
	if(.not.kas_vaikselt) print*, "All done:D", t2-t1

end program gm