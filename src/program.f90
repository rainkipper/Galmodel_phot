program gm
	use fill_comp_image_module
	use	read_input_module
	use fitting_multinest_module
	use only_output_image_module

	type(all_comp_type) :: all_comp, all_comp_input
	type(comp_input_type), dimension(:), allocatable :: input_comp
	type(comp_image_real_type) :: mdl1, mdl2, mdl
	type(image_type), dimension(:), allocatable :: images
	type(filter_type), dimension(:), allocatable :: filters
	integer, parameter :: N_comp = 1
	integer :: test_image_number
	character(len=*), parameter :: fname="comp_im.txt"
	character(len=default_character_length), parameter :: input_comp_file = "Input/input_comp.txt"
	character(len=default_character_length), parameter :: input_image_file = "Input/input_images.txt"
	logical :: kas_comp_im, kas_los
	integer :: i, N
	real(rk), dimension(:,:), allocatable :: pilt
	real(rk) :: t1,t2 !aja mootmine
	real(rk) :: ll
	logical, parameter :: kas_fitib_vs_lihtsalt_ouput = .not..false.
	
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
		call asenda_viited(input_comp, all_comp) !siin all on ka init comp
		call create_output_image(all_comp, images(1), pilt)
		call write_matrix_to_fits(pilt, images(1)%output_mdl_file)
	end if
	

! 	ll = calc_log_likelihood(all_comp, images)
	
	call  cpu_time(t2)
	print*, "All done:D", t2-t1

contains
    
end program gm