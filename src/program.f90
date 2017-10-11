program gm
	use fill_comp_image_module
	use	read_input_module
	use fitting_module
	use only_output_image_module
	use mynr_amoeba_module

	type(all_comp_type) :: all_comp !, all_comp_input
	type(comp_input_type), dimension(:), allocatable :: input_comp
! 	type(comp_image_real_type) :: mdl1, mdl2, mdl
	type(image_type), dimension(:), allocatable :: images
	type(filter_type), dimension(:), allocatable :: filters
	integer, parameter :: N_comp = 1
! 	integer :: i, iter
! 	integer :: test_image_number
	character(len=*), parameter :: fname="comp_im.txt"
! 	character(len=default_character_length), parameter :: input_comp_file = "Input_sdss/input_comp.txt"
! 	character(len=default_character_length), parameter :: input_image_file = "Input_sdss/input_images.txt"
	character(len=default_character_length), parameter :: input_comp_file = "Input/Mock/input_comp.txt"
	character(len=default_character_length), parameter :: input_param_file = "Input/Mock/input_parameters.txt"
	character(len=default_character_length), parameter :: input_image_file = "Input/Mock/input_images.txt"
! 	logical :: kas_comp_im, kas_los
! 	integer :: i, N
	real(rk), dimension(:,:), allocatable :: pilt
! 	real(rk), dimension(:), allocatable :: vec, vec2
	real(rk) :: t1,t2 !aja mootmine
! 	real(rk) :: ll
	logical, parameter :: kas_fitib_vs_lihtsalt_ouput = .not..false.
	call random_seed()
	call cpu_time(t1)
	call read_constants_and_parameters(input_param_file)
	
	!
! 	allocate(pilt(1:10, 1:9))
! 	allocate(vec(1:10))
! 	allocate(vec2(1:9))
! 	call random_number(pilt); pilt = pilt * 30
! 	do i=1,10;
! 		vec2 = pilt(i,:)
! 		vec(i)=v2hendatav(vec2) ;
! 	end do
! 	print*, vec
! 	call my_amoeba(pilt, vec, 0.001_rk, v2hendatav, iter)
! 	print*, "iter", iter
! 	print*, vec
! 	print*, pilt(1,:)
! 	stop "amoeba test tehtud"
	
	!
	! ================== mingi testgalaktika tekitamine =================
	!
	call read_components(input_comp_file, input_comp)
	call read_images_and_filters(input_image_file, images, filters)
	call init_images(images) !peamiselt psf crop

	if(kas_fitib_vs_lihtsalt_ouput) then
		call fit_galaxy(images, input_comp, all_comp)
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
!     function v2hendatav(x) result(res)
! 		!amoobi testimiseks
! 		implicit none
! 		real(rk), intent(in), dimension(:) :: x
! 		integer :: i
! 		real(rk) :: res
! 		res = 0.0
! 		do i=1,size(x,1)
! 			res = res + (x(i)-i)**2
! 		end do
! 	end function v2hendatav
end program gm