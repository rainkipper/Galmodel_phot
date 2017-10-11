module fitting_module
	use fitting_multinest_module
	use fitting_amoeba_module

	
contains
		subroutine fit_galaxy(images, input_comps, all_comp)
			implicit none
			type(comp_input_type), dimension(:), allocatable, intent(inout), target :: input_comps
			type(all_comp_type), intent(out) :: all_comp !ehk v2ljundiks
			type(image_type), dimension(:), allocatable, intent(in) :: images
			
			!
			! initsialiseerimise asjad, et teha all_comp jm..
			! sisuliselt, yhine koikidel fittijatel
			!
			call convert_input_comp_to_all_comp(input_comps, all_comp)
			call asenda_viited(input_comps, all_comp)
			all_comp%comp(:)%adaptive_image_number = -1 !default -1, et teeks uue adaptiivse pildi esimene kord
			call init_calc_log_likelihood(all_comp, images) !s2ttib likelihoodi mooduli muutujad, et v2hendada arvutamisis
			
			

! 			call fittimine_multinest(images, input_comps, all_comp)
			call fittimine_amoeba(images, input_comps, all_comp)
			
			
		end subroutine fit_galaxy
end module fitting_module