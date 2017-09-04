module filters_module
	use constants_module
	type filter_type
		character(len=default_character_length) 							:: name
		character(len=default_character_length), dimension(:), allocatable 	:: population_names
		real(rk), dimension(:), allocatable 								:: population_mass_to_light_ratios
		real(rk) 															:: Mag_sun !abs mag
		real(rk)															:: ZP 	
		contains
			procedure :: calc_ML_ratio
	end type filter_type
contains
	function calc_ML_ratio(filter, dist, population_name) result(res)
		implicit none
		class(filter_type), intent(in) :: filter
		character(len=default_character_length), intent(in) :: population_name
		real(rk), intent(in) :: dist
		real(rk) :: res
		real(rk) :: pop_ML
		integer :: i
		pop_ML = -1234.5
		do i=1,size(filter%population_names, 1)
			if(trim(population_name) == trim(filter%population_names(i))) then
				pop_ML = filter%population_mass_to_light_ratios(i)
			end if
		end do
		if(pop_ML == -1234.5) then
			print*, "Not possible to match population in filter section"
			stop
		end if
		res = 10**(0.4*(filter%ZP-filter%Mag_sun) + 6.0 - 2.0*log10( dist ) ) / pop_ML
	end function  calc_ML_ratio
	subroutine create_test_filters(filters)
		implicit none
		type(filter_type), dimension(:), allocatable :: filters
		integer :: N
		N = 2
		allocate(filters(1:N))
		
		filters(1)%name = "HST_v"
		filters(1)%population_names = ["young_population", "old_population"]
		filters(1)%population_mass_to_light_ratios = [1.0, 0.1]
		filters(1)%Mag_sun = 4.74_rk
		filters(1)%ZP = 26.50512_rk
		
		filters(2)%name = "HST_z"
		filters(2)%population_names = ["young_population", "old_population"]
		filters(2)%population_mass_to_light_ratios = [0.1, 1.0]
		filters(2)%Mag_sun = 5.94_rk
		filters(2)%ZP = 24.86663_rk
	end subroutine create_test_filters
end module filters_module