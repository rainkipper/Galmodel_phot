module los_real_integration_module
	use integration_module
	use constants_module
	use comp_module
	real(rk) :: kauguse_piir = 20.0_rk !integreerimise rada
contains
	!
	!keerulisemate systeemide korral (kus mitut objekti korraga fititakse), on vaja koike teha fyysikalistes koordinaatides, mitte komponendi koordinaatides ning siis integreerimine_los peab olema fyysikalises ning ruum_ptr peaks sisendi saama fyysikalistes koordinaatides. 
	!
	function integreerimine_los(ruum_ptr, comp, Xc, Yc) result(res)
		implicit none
		interface
			function ruum_ptr(R,z,theta) result(res)
				import rk
				real(rk), intent(in) :: R,z
				real(rk), intent(in), optional :: theta
				real(rk) :: res
			end function ruum_ptr
		end interface
		type(comp_type), intent(in) :: comp
		real(rk), intent(in) :: Xc, Yc
		real(rk) :: res
		real(rk) :: integreerimise_t2psus !

		integreerimise_t2psus = ruum_ptr(Xc, Yc)*1.0e-3
		res = integrate(los, -1.0*kauguse_piir, kauguse_piir, integreerimise_t2psus)
	contains
		function los(d) result(res)
			implicit none
			real(rk), intent(in) :: d
			real(rk) :: res
			real(rk) :: R, z, theta
			call XcYcl_to_Rztheta(Xc, Yc, d,  comp%sin_incl, comp%cos_incl, comp%tan_incl, comp%sec_incl, comp%theta0, R, z, theta)
			res = ruum_ptr(R,z,theta)
		end function los
	end function integreerimine_los
	
! 	function scalar_over_los(func, kauguse_piir) result(res)
! 		implicit none
! 		interface
! 			function func(d) result(res)
! 				import rk
! 				real(rk), intent(in) :: d
! 				real(rk) :: res
! 			end function func
! 		end interface
! 		real(rk), intent(in) :: kauguse_piir
! 		!type(accuracy_type) :: t2psus
! 		real(rk) :: t2psus
! 		real(rk) :: res
! 		real(rk), allocatable, dimension(:,:) :: fvalue
! 		integer :: errcode
! 		t2psus = 1.0e-4*func(0.0_rk)
! ! 		integrate_gauss_adaptive(func,a,b,acc,nmax,nmin,errcode,incracc,fvalue,log_base)
! ! 		res = integrate_gauss_adaptive(func=func, a=-1.0*kauguse_piir, b=kauguse_piir, acc=t2psus, nmax=10, nmin=2, errcode=errcode, incracc = .true.)
! 		res = integrate(func, -1.0*kauguse_piir, kauguse_piir, t2psus)
! 	end function scalar_over_los
end module los_real_integration_module