module mass_to_lum_image_module
	use comp_image_module
	use comp_module
	use populations_module
	use dust_module
	use filters_module
	private
	public :: make_lum_image, init_lum_dependencies
	interface make_lum_image
		module procedure :: ML_abil_image
		module procedure :: population_dust_image
	end interface make_lum_image
	type(distribution_1D_type) :: sol_spec !
contains
	subroutine init_lum_dependencies(sol_file_name, populations_file_name)
		implicit none
		character(len = default_character_length), intent(in) :: sol_file_name, populations_file_name
		call read_solar_spectrum()
		call create_sdss_filters() !see peaks muutuma mingi aeg
		call init_populations(populations_file_name)
	contains
		subroutine read_solar_spectrum()
			implicit none
			real(rk), dimension(:,:), allocatable :: toores
			integer :: N
			toores = read_tabel(sol_file_name, 4)
			toores(:,2) = toores(:,2) * 100.0 * tan(1.0/3600.0*pi/180.0/10.0_rk)**2
			!rida 3 on tilt, mille kohta ei tea midagi
			toores(:,4) = toores(:,4) * 100.0 * tan(1.0/3600.0*pi/180.0/10.0_rk)**2
			toores(:,1) = 10.0*toores(:,1) !ehk koik on angstromid
			allocate(sol_spec%x(1:size(toores, 1)))
			allocate(sol_spec%f(1:size(toores, 1)))
			sol_spec%x = toores(:,1)
			sol_spec%f = toores(:,2)
		end subroutine read_solar_spectrum
		
	end subroutine init_lum_dependencies
	
	function ML_abil_image(comp_im, filter, dist) result(res)
		implicit none
		type(comp_image_real_type), intent(in) :: comp_im
		type(filter_type), intent(in) :: filter
		real(rk), intent(in) :: dist
		real(rk) :: ML
		real(rk), dimension(:,:), allocatable :: res
		if(.not.allocated(res)) allocate(res(1:size(comp_im%mx,1), 1:size(comp_im%mx,2)))
		res = 0.0
	end function ML_abil_image
	
	function population_dust_image(comp_im, filter, pop, dust, dist, comp_nimi) result(res)
		implicit none
		type(comp_image_real_type), intent(in) :: comp_im
		class(comp_type), intent(in) :: dust
		type(filter_type), intent(in) :: filter
		type(population_type), intent(in) :: pop
		character(len=default_character_length), intent(in), optional :: comp_nimi
		real(rk) :: tau0, tau
		real(rk), intent(in) :: dist
		real(rk) :: ML
		real(rk), dimension(:), allocatable :: lainepikkus
		real(rk) ::  sollum, dustlesslum, dustylum
		real(rk), dimension(:,:), allocatable :: res
		integer :: Nlambda
		integer :: i
		real(rk) :: minraja, maxraja
		real(rk) :: ML_tolmuga, ML_tolmuta
		
		if(.not.allocated(res)) allocate(res(1:size(comp_im%mx,1), 1:size(comp_im%mx,2)))
		res = 0.0
		!
		! integreerimised ML jaoks yle lainepikkuste
		!
		minraja = max(minval(sol_spec%x, 1), minval(pop%spec%x, 1), minval(filter%spec%x, 1))
		maxraja = min(maxval(sol_spec%x, 1), maxval(pop%spec%x, 1), maxval(filter%spec%x, 1))
		Nlambda = (maxraja-minraja)/lainepikkuse_dl !dl on voetud const juurest
		allocate(lainepikkus(1:Nlambda))
		lainepikkus(1) = minraja
		do i=2,Nlambda
			lainepikkus(i) = lainepikkus(i-1)+lainepikkuse_dl
		end do
		call dust%prof_den%get_val("tau0", tau0) !kordaja
		tau = tau0 / dust%cos_incl !kaldenurga parand
! 		print*, "CP0"
! 		print*, "min ja max"
! 		print*, minval(sol_spec%x, 1), minval(pop%spec%x, 1), minval(filter%spec%x, 1)
! 		print*, maxval(sol_spec%x, 1), maxval(pop%spec%x, 1), maxval(filter%spec%x, 1)
		sollum = lainepikkuse_dl * sum(interpolate_1D(sol_spec, lainepikkus) * interpolate_1D(filter%spec, lainepikkus)) !tavaline Riemanni integraal
		dustlesslum = lainepikkuse_dl * sum(interpolate_1D(pop%spec, lainepikkus) * interpolate_1D(filter%spec, lainepikkus))
		!kui tahta mingit muud tolmu varianti panna (sh 2mikroni feature-ga), siis see tuleks panna j2rgneva rea peale exp alla.
		dustylum = lainepikkuse_dl * sum((interpolate_1D(pop%spec, lainepikkus)) * interpolate_1D(filter%spec, lainepikkus) * exp(-abs(tau)*dust_kappa(lainepikkus))) 
! 		print*, "sol spec", sol_spec%x(20), sol_spec%f(20)
! 		print*, "sol na", count(isnan(interpolate_1D(sol_spec, lainepikkus))), sum((interpolate_1D(sol_spec, lainepikkus)))
! 		print*, "tau na", count(isnan(exp(-abs(tau0)*dust_kappa(lainepikkus)))), sum(exp(-abs(tau0)*dust_kappa(lainepikkus)))
! 		print*, "filt na", count(isnan(interpolate_1D(filter%spec, lainepikkus))), sum((interpolate_1D(filter%spec, lainepikkus)))
! 		print*, "pop na", count(isnan(interpolate_1D(pop%spec, lainepikkus))), sum((interpolate_1D(pop%spec, lainepikkus)))
! 		print*, "lambda", sum(lainepikkus)
! 		print*, "lainepikkuse_dl",lainepikkuse_dl, Nlambda
! 		print*, "sol spec piirid", minval(sol_spec%x, 1), maxval(sol_spec%x, 1)
! 		print*, "rajad", minraja, maxraja
! 		print "(A,A,A,2F10.5, E15.8)", trim(filter%name), " ", trim(pop%name), sollum/dustlesslum, sollum/dustylum, sollum
		!
		! ML arvutamine
		!

! 		print*, "tau = ", tau, tau0
		ML_tolmuta = 1.0/(dustlesslum/sollum) !mass on 1
		ML_tolmuga = 1.0/(dustylum/sollum)
! 		print*, "ML_tolmu(g/t)a", ML_tolmuga, ML_tolmuta
		!
		! pildi peale rakendamine
		!
		if(tau0<0) then
			res = comp_im%M_enne_tasandit/ML_tolmuga + comp_im%M_p2rast_tasandit/ML_tolmuta !taul pole miinust poole definitsiooni p2rast
		else
			res = comp_im%M_enne_tasandit/ML_tolmuta + comp_im%M_p2rast_tasandit/ML_tolmuga
		end if
		res = res * filter%mass_to_obs_unit(dist) !pildi yhikuteks... seda voib ka varem teisendada
		if(count(isnan(res))>0) then
			print "(A,A,A,A)", "Tekkisid mingid NA-d (filt, pop):", trim(filter%name), " ", trim(pop%name)
		end if
! 		if(present(comp_nimi)) then
! 			call write_matrix_to_fits(res, "Output/comp_im_"//trim(comp_nimi)//"_"//trim(filter%name)//".fits")
! 		end if
		
! 		print*, "ML", ML_tolmuta, ML_tolmuga
! 		print*, "CP inf"
	end function population_dust_image
	
	
end module mass_to_lum_image_module