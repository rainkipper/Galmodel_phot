module only_output_image_module
	use images_module
	use fill_comp_image_module
	use likelihood_module !vajalik massi t2psuse saamiseks
! 	logical, parameter, private :: kas_los = .true.
	logical, parameter, private :: kas_comp_im = .true.
contains
	subroutine create_output_image(all_comp, image, pilt)
		type(all_comp_type), intent(inout) :: all_comp
		type(image_type), intent(in) :: image
		type(image_type), dimension(:), allocatable :: im
		type(comp_image_real_type), dimension(:), allocatable :: mudelid
		real(rk), dimension(:), allocatable :: weights
		real(rk), dimension(:,:), allocatable, intent(out) :: pilt
		integer :: i
! 		!vajaliku abs t2psuse leidmine... optional
! 		allocate(im(1:1)); im(1) = image !jama, et saada funktsiooni sisendiks oiget kuju
! 		do i=1,all_comp%N_comp
! 			all_comp%comp(i)%mass_abs_tol = leia_massi_abs_tol(all_comp%comp(i), im)
! ! 			all_comp%comp(i)%mass_abs_tol = 1.0e10 !ehk ei arvesta
! 		end do
! 		!mudelpildi tegemine
! 		allocate(weights(1:all_comp%N_comp))
! 		allocate(mudelid(1:all_comp%N_comp))
! 		do i=1,size(mudelid, 1)
! 			call create_comp_image_from_obs(mudelid(i), image)
! 			call fill_comp_image(all_comp, i, mudelid(i), kas_comp_im, kas_los)
! 			weights(i) = image%filter%mass_to_obs_unit(all_comp%comp(i)%dist, all_comp%comp(i)%population_name)
! 		end do
! 		pilt = combine_comp_images_to_make_image(mudelid, weights)
	end subroutine create_output_image
end module only_output_image_module