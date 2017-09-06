module filters_module
	use constants_module

	type filter_type
		character(len=default_character_length) 							:: name
		character(len=default_character_length), dimension(:), allocatable 	:: population_names
		real(rk), dimension(:), allocatable 								:: population_mass_to_light_ratios
		real(rk) 															:: Mag_sun !abs mag
		real(rk)															:: ZP 	
		procedure(teisenduse_kuju), pointer :: mass_to_obs_unit
	end type filter_type
	interface
		function teisenduse_kuju(filter, dist, population_name)
			import rk
			import default_character_length
			import filter_type
			real(rk) :: teisenduse_kuju
			real(rk), intent(in) :: dist
			class(filter_type), intent(in) :: filter
			character(len=default_character_length), intent(in) :: population_name
		end function teisenduse_kuju
	end interface
	interface create_test_filters
! 		module procedure :: create_hst_filters
		module procedure :: create_sdss_filters
	end interface create_test_filters

contains
	function calc_counts_mass_ratio(filter, dist, population_name) result(res)
		!tulemusega korrutades saab massist countsid
		implicit none
		class(filter_type), intent(in) :: filter
		character(len=default_character_length), intent(in) :: population_name
		real(rk), intent(in) :: dist
		real(rk) :: res
		real(rk) :: pop_ML
		integer :: i

		!M/L osa, mis tuleneb otse populatsioonist
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
		!pildi enda yhikute arvestamine (counts eeldusel)
		res = 10**(0.4*(filter%ZP-filter%Mag_sun) + 6.0 - 2.0*log10( dist ) ) / pop_ML
	end function  calc_counts_mass_ratio
	subroutine create_hst_filters(filters)
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
		filters(1)%mass_to_obs_unit => calc_counts_mass_ratio
		
		filters(2)%name = "HST_z"
		filters(2)%population_names = ["young_population", "old_population"]
		filters(2)%population_mass_to_light_ratios = [0.1, 1.0]
		filters(2)%Mag_sun = 5.94_rk
		filters(2)%ZP = 24.86663_rk
		filters(2)%mass_to_obs_unit => calc_counts_mass_ratio
	end subroutine create_hst_filters
	
	
		
	subroutine create_sdss_filters(filters)
! 		    filter	Vega	AB	  Apparent_AB mag_Vega
! 		26	SDSS u'	5.46	6.45	-26.12	0.995
! 		27	SDSS g'	5.22	5.14	-26.36	-0.087
! 		28	SDSS r'	4.50	4.65	-27.08	0.147
! 		29	SDSS i'	4.16	4.54	-27.42	0.376
! 		30	SDSS z'	4.01	4.52	-27.58	0.518
		implicit none
		type(filter_type), dimension(:), allocatable :: filters
		integer :: N
		N = 5
		allocate(filters(1:N))
		
		filters(1)%name = "sdss_u"
		filters(1)%population_names = ["young_population", "old_population"]
		filters(1)%population_mass_to_light_ratios = [1.0, 0.1]
		filters(1)%Mag_sun = 6.45_rk
		filters(1)%ZP = 22.5_rk !yhikute ning ABmag erip2ra
		filters(1)%mass_to_obs_unit => calc_counts_mass_ratio
		
		
		

	end subroutine create_sdss_filters
end module filters_module