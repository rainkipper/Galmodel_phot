module fitting_module
	use fitting_multinest_module
	use fitting_amoeba_module
	use fitting_omalooming_module
	!iga fittija peab saama sisse pildid, input_comps ning v2ljundis all_comp (kuigi seda otseselt vaja pole)
	!iga fittija peab lopus tegema parima  input_comps-i, mis l2heb outputi tegemiseks
	
contains
		subroutine fit_galaxy(images, input_comps, all_comp)
			implicit none
			type(comp_input_type), dimension(:), allocatable, intent(inout), target :: input_comps
			type(all_comp_type), intent(out) :: all_comp !ehk v2ljundiks
			type(image_type), dimension(:), allocatable, intent(in) :: images
			real(rk), dimension(:,:), allocatable :: punktid
			real(rk) :: LL
			!
			! initsialiseerimise asjad, et teha all_comp jm..
			! sisuliselt, yhine koikidel fittijatel
			!
			call convert_input_comp_to_all_comp(input_comps, all_comp)
			call asenda_viited(input_comps, all_comp)
			all_comp%comp(:)%adaptive_image_number = -1 !default -1, et teeks uue adaptiivse pildi esimene kord
			call init_calc_log_likelihood(all_comp, images) !s2ttib likelihoodi mooduli muutujad, et v2hendada arvutamisis
			select case( mis_fittimise_algoritm)
			case(0)
				if(.not.kas_vaikselt) print*, "Only output"
				call convert_input_comp_to_all_comp(input_comps, all_comp)
				call asenda_viited(input_comps, all_comp) !all_comp muutujas asendamine
				LL =  calc_log_likelihood(all_comp, images)
			case(1)
				if(.not.kas_vaikselt) print*, "Multinest fittimine"
				call fittimine_multinest(images, input_comps, all_comp)
			case(2)
				if(.not.kas_vaikselt) print*, "Amoeba fittimine"
				call fittimine_amoeba(images, input_comps, all_comp)
			case(3)
				if(.not.kas_vaikselt) print*, "Multinest+amoeba fittimine"
				call fittimine_multinest(images, input_comps, all_comp)
				call read_points_for_amoeba(multinest_output_header, punktid)
				call fittimine_amoeba(images, input_comps, all_comp, punktid)
			case(4)
				if(.not.kas_vaikselt) print*, "Amoeba from previous multinest"
				call read_points_for_amoeba(multinest_output_header, punktid)
				call fittimine_amoeba(images, input_comps, all_comp, punktid)
			case(5)
				if(.not.kas_vaikselt) print*, "fittimine omaloominguga"
				call fittimine_omalooming(images, input_comps, all_comp)
			case default
				stop "Niisugust fittimise meetodit pole olemas"
			end select
			
			if(.not.kas_vaikselt) print*, "Kokku oli LL arvutusi", LL_counter	
			call output_images(input_comps, images)
			call output_like_input(input_comps)
			if(mis_fittimise_tyyp == 2) call output_ML(input_comps, images)
			
		contains
			subroutine read_points_for_amoeba(file_path, res)
				implicit none
				character(len=*), intent(in) :: file_path
				real(rk), dimension(:,:), allocatable, intent(inout) :: res
				integer :: iunit, N, i
				N = leia_vabade_parameetrite_arv(input_comps) !moodulist comp.f90
				if(allocated(res)) deallocate(res); allocate(res(1:(N+1),1:N)); res = -1.0 !tulemuse tekitamine
				iunit = 21
				open(file=file_path//"phys_live.points", unit=iunit, action = "read")
				do i = 1,N+1
					read(fmt=*, unit=iunit) res(i,:)
				end do
				close(unit=iunit)
			end subroutine read_points_for_amoeba
		end subroutine fit_galaxy
end module fitting_module