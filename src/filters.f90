module filters_module
	use constants_module
	use distributions_module
	use file_operations_module
	type filter_type
		character(len=default_character_length) 							:: name
		character(len=default_character_length), dimension(:), allocatable 	:: population_names
		real(rk), dimension(:), allocatable 								:: population_mass_to_light_ratios
		real(rk), dimension(:), allocatable 								:: population_tau_kordajad
		real(rk) 															:: Mag_sun !abs mag
		real(rk)															:: ZP 	
		type(distribution_1D_type)												:: spec !vaja ML suhete arvutamisel
		procedure(filtri_yhikute_teisenduse_kuju), pointer :: mass_to_obs_unit
! 		procedure(tau_kordaja_saamise_kuju), pointer :: get_tau_kordaja
		
	end type filter_type
	interface
		function filtri_yhikute_teisenduse_kuju(filter, dist, population_name) result(res)
			import rk
			import default_character_length
			import filter_type
			real(rk) :: res
			real(rk), intent(in) :: dist
			class(filter_type), intent(in) :: filter
			character(len=default_character_length), intent(in), optional :: population_name
		end function filtri_yhikute_teisenduse_kuju
		function tau_kordaja_saamise_kuju(filter, population_name) result(res)
			import rk
			import default_character_length
			import filter_type
			real(rk) :: res
			class(filter_type), intent(in) :: filter
			character(len=default_character_length), intent(in), optional :: population_name
		end function tau_kordaja_saamise_kuju
	end interface
	interface create_test_filters
! 		module procedure :: create_hst_filters
		module procedure create_sdss_filters
! 			subroutine create_sdss_filters(filters)
! 				import filter_type
! 				class(filter_type), dimension(:), allocatable, intent(out) :: filters
! 			end subroutine create_sdss_filters
	end interface
	type(filter_type), dimension(:), allocatable, target :: filters
contains
	
	
	
! 	function get_tau_kordaja(filter, population_name) result(res)
! 		!tulemusega korrutades saab massist countsid
! 		implicit none
! 		class(filter_type), intent(in) :: filter
! 		character(len=default_character_length), intent(in), optional :: population_name
! 		real(rk) :: res
! 		integer :: i
!
! 		!M/L osa, mis tuleneb otse populatsioonist
! 		if(present(population_name)) then
! 			res = -1234.5
! 			do i=1,size(filter%population_names, 1)
! 				if(trim(population_name) == trim(filter%population_names(i))) then
! 					res = filter%population_tau_kordajad(i)
! 				end if
! 			end do
! 			if(res == -1234.5) then
! 				print*, "Not possible to match population in filter section (tau)"
! 				stop
! 			end if
! 		end if
! 	end function  get_tau_kordaja

	
function calc_counts_mass_ratio(filter, dist, population_name) result(res)
	!tulemusega korrutades saab massist countsid
	implicit none
	class(filter_type), intent(in) :: filter
	character(len=default_character_length), intent(in), optional :: population_name
	real(rk), intent(in) :: dist
	real(rk) :: res
	real(rk) :: pop_ML
	integer :: i

	!M/L osa, mis tuleneb otse populatsioonist
	if(present(population_name)) then
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
	else
		pop_ML = 1.0
	end if
	!pildi enda yhikute arvestamine (counts eeldusel)
	res = 10**(0.4*(filter%ZP-filter%Mag_sun) + 6.0 - 2.0*log10( dist ) ) / pop_ML
end function  calc_counts_mass_ratio


	
	subroutine create_hst_filters()
		implicit none
! 		type(filter_type), dimension(:), allocatable :: filters
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
	
	
		
	subroutine create_sdss_filters_test()
		!Mass heledus suhted on arvutatud mingi SSP pohjal Bruzual-Charlot 
		!http://www.bruzual.org/bc03/Updated_version_2016/
		!R script ML_arvutamine.R t2psete numbrite saamiseks
! 		    filter	Vega	AB	  Apparent_AB mag_Vega
! 		26	SDSS u'	5.46	6.45	-26.12	0.995
! 		27	SDSS g'	5.22	5.14	-26.36	-0.087
! 		28	SDSS r'	4.50	4.65	-27.08	0.147
! 		29	SDSS i'	4.16	4.54	-27.42	0.376
! 		30	SDSS z'	4.01	4.52	-27.58	0.518
		implicit none
! 		type(filter_type), dimension(:), allocatable :: filters
		integer :: N
		N = 5
		allocate(filters(1:N))
		
		filters(1)%name = "SDSS_u"
		filters(1)%population_names = ["pop_logt6", "pop_logt9", "pop_logt10"]
		filters(1)%population_mass_to_light_ratios = [0.1591402, 1.5185881, 3.498968]
		filters(1)%Mag_sun = 6.45_rk
		filters(1)%ZP = 22.5_rk !yhikute ning ABmag erip2ra
		filters(1)%mass_to_obs_unit => calc_counts_mass_ratio
		
		
		filters(2)%name = "SDSS_g"
		filters(2)%population_names = ["pop_logt6", "pop_logt9", "pop_logt10"]
		filters(2)%population_mass_to_light_ratios = [0.2861805, 1.4190575, 3.546990]
		filters(2)%Mag_sun = 5.14_rk
		filters(2)%ZP = 22.5_rk !yhikute ning ABmag erip2ra
		filters(2)%mass_to_obs_unit => calc_counts_mass_ratio
		
		filters(3)%name = "SDSS_r"
		filters(3)%population_names = ["pop_logt6", "pop_logt9", "pop_logt10"]
		filters(3)%population_mass_to_light_ratios = [0.4305861, 1.3543783, 3.312836]
		filters(3)%Mag_sun = 4.65_rk
		filters(3)%ZP = 22.5_rk !yhikute ning ABmag erip2ra
		filters(3)%mass_to_obs_unit => calc_counts_mass_ratio
		
		filters(4)%name = "SDSS_i"
		filters(4)%population_names = ["pop_logt6", "pop_logt9", "pop_logt10"]
		filters(4)%population_mass_to_light_ratios = [0.5156797, 1.2004387, 2.941034]
		filters(4)%Mag_sun = 4.54_rk
		filters(4)%ZP = 22.5_rk !yhikute ning ABmag erip2ra
		filters(4)%mass_to_obs_unit => calc_counts_mass_ratio
		
		filters(5)%name = "SDSS_z"
		filters(5)%population_names = ["pop_logt6", "pop_logt9", "pop_logt10"]
		filters(5)%population_mass_to_light_ratios = [0.5834451, 0.9559606, 2.369627]
		filters(5)%Mag_sun = 4.52_rk
		filters(5)%ZP = 22.5_rk !yhikute ning ABmag erip2ra
		filters(5)%mass_to_obs_unit => calc_counts_mass_ratio
		
		
		

	end subroutine create_sdss_filters_test
	
	subroutine create_sdss_filters()
		!Blantoni templateite pohjal, kus on tolm ja spektri jooned ka olemas
! 		    filter	Vega	AB	  Apparent_AB mag_Vega
! 		26	SDSS u'	5.46	6.45	-26.12	0.995
! 		27	SDSS g'	5.22	5.14	-26.36	-0.087
! 		28	SDSS r'	4.50	4.65	-27.08	0.147
! 		29	SDSS i'	4.16	4.54	-27.42	0.376
! 		30	SDSS z'	4.01	4.52	-27.58	0.518
		implicit none
		character(len=default_character_length) :: path, file
		real(rk), dimension(:,:), allocatable :: mx
! 		type(filter_type), dimension(:), allocatable :: filters
		character(len=default_character_length), dimension(1:5) :: filenames = ["u.dat.txt", "g.dat.txt", "r.dat.txt", "i.dat.txt", "z.dat.txt"]
		integer :: N, i
		N = 5
		allocate(filters(1:N))

		filters(1)%population_names = ["Blanton_1", "Blanton_2", "Blanton_3", "Blanton_4", "Blanton_5"]
		filters(2)%population_names = ["Blanton_1", "Blanton_2", "Blanton_3", "Blanton_4", "Blanton_5"]
		filters(3)%population_names = ["Blanton_1", "Blanton_2", "Blanton_3", "Blanton_4", "Blanton_5"]
		filters(4)%population_names = ["Blanton_1", "Blanton_2", "Blanton_3", "Blanton_4", "Blanton_5"]
		filters(5)%population_names = ["Blanton_1", "Blanton_2", "Blanton_3", "Blanton_4", "Blanton_5"]
		filters(1)%ZP = 22.5_rk !yhikute ning ABmag erip2ra
		filters(2)%ZP = 22.5_rk !yhikute ning ABmag erip2ra
		filters(3)%ZP = 22.5_rk !yhikute ning ABmag erip2ra
		filters(4)%ZP = 22.5_rk !yhikute ning ABmag erip2ra
		filters(5)%ZP = 22.5_rk !yhikute ning ABmag erip2ra
		
		filters(1)%name = "SDSS_u"
		filters(1)%population_mass_to_light_ratios = [23.3301748622864,0.00121417016930545,0.439758219586101,7.19578436927323,0.148393334920314]
		filters(1)%population_tau_kordajad = [5.0, 4.0, 3.0, 2.0, 1.0] !testiks
		filters(1)%Mag_sun = 6.45_rk
		filters(1)%mass_to_obs_unit => calc_counts_mass_ratio
! 		filters(1)%get_tau_kordaja=>get_tau_kordaja
		
		
		filters(2)%name = "SDSS_g"
		filters(2)%population_mass_to_light_ratios = [8.22717332271278,0.00432318253378151,0.454359561061207,5.00051509903871,0.255755562862128]
		filters(2)%population_tau_kordajad = [5.0, 4.0, 3.0, 2.0, 1.0] !testiks
		filters(2)%Mag_sun = 5.14_rk
		filters(2)%mass_to_obs_unit => calc_counts_mass_ratio
! 		filters(2)%get_tau_kordaja=>get_tau_kordaja
		
		filters(3)%name = "SDSS_r"
		filters(3)%population_mass_to_light_ratios = [3.95596422691167,0.00883400910870545,0.495677538251191,3.85409345140553,0.339372528787265]
		filters(3)%population_tau_kordajad = [5.0, 4.0, 3.0, 2.0, 1.0] !testiks
		filters(3)%Mag_sun = 4.65_rk
		filters(3)%mass_to_obs_unit => calc_counts_mass_ratio
! 		filters(3)%get_tau_kordaja=>get_tau_kordaja
		
		filters(4)%name = "SDSS_i"
		filters(4)%population_mass_to_light_ratios = [2.44997274893593,0.0158117700010423,0.562029102838556,3.12625952921625,0.344491806973725]
		filters(4)%population_tau_kordajad = [5.0, 4.0, 3.0, 2.0, 1.0] !testiks
		filters(4)%Mag_sun = 4.54_rk
		filters(4)%mass_to_obs_unit => calc_counts_mass_ratio
! 		filters(4)%get_tau_kordaja=>get_tau_kordaja
		
		filters(5)%name = "SDSS_z"
		filters(5)%population_mass_to_light_ratios = [1.54877795627825,0.0184719568573701,0.512100014012443,2.51197061658436,0.357531160626106]
		filters(5)%population_tau_kordajad = [5.0, 4.0, 3.0, 2.0, 1.0] !testiks
		filters(5)%Mag_sun = 4.52_rk
		filters(5)%mass_to_obs_unit => calc_counts_mass_ratio
! 		filters(5)%get_tau_kordaja=>get_tau_kordaja

		path = "/Users/rain/Suured_tegemised/Califa_SDSS/Populations/"
		do i=1,N
			file = repeat(" ", default_character_length); file = trim(path)//trim(filenames(i))
			if(allocated(mx)) deallocate(mx); mx = read_tabel(file, 5)
			allocate(filters(i)%spec%x(1:size(mx, 1)))
			allocate(filters(i)%spec%f(1:size(mx, 1)))
			filters(i)%spec%x = mx(:,1)
			filters(i)%spec%f = mx(:,4)
		end do
		



	end subroutine create_sdss_filters
end module filters_module