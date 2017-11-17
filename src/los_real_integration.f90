module los_real_integration_module
	use integration_module
	use comp_module
contains
	!
	!keerulisemate systeemide korral (kus mitut objekti korraga fititakse), on vaja koike teha fyysikalistes koordinaatides, mitte komponendi koordinaatides ning siis integreerimine_los peab olema fyysikalises ning ruum_ptr peaks sisendi saama fyysikalistes koordinaatides. 
	!
	function integreerimine_los(ruum_ptr, comp, Xc, Yc, lopmatus) result(res)
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
		real(rk), intent(in) :: Xc, Yc, lopmatus
		real(rk) :: res
		real(rk) :: integreerimise_t2psus !

		integreerimise_t2psus = comp%massi_abs_tol_los
! 		print*, "los t2psus", integreerimise_t2psus
		res = integrate(func = los, a = -1.0*lopmatus, b = lopmatus, acc = integreerimise_t2psus, nmin = integration_min_iter, nmax = integration_max_iter)
	contains
		function los(d) result(res)
			implicit none
			real(rk), intent(in) :: d
			real(rk) :: res
			real(rk) :: R, z, theta
			call XcYcl_to_Rztheta(Xc, Yc, d,  comp%sin_incl, comp%cos_incl, comp%tan_incl, comp%sec_incl, comp%theta0, R, z, theta)
			res = ruum_ptr(R,z,theta)
! 			print*, d,res
		end function los
	end function integreerimine_los
	
	function integreerimine_los_dustplane(ruum_ptr, comp, Xc, Yc, lopmatus, kas_enne_tasandit) result(res)
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
		real(rk), intent(in) :: Xc, Yc, lopmatus
		logical, intent(in) :: kas_enne_tasandit
		real(rk) :: res
		real(rk) :: integreerimise_t2psus !

		integreerimise_t2psus = comp%massi_abs_tol_los
! 		print*, "los t2psus", integreerimise_t2psus
		!res1 ja res2 on heledused vastavalt enne ja p2rast tolmu tasandit
		if(kas_enne_tasandit) then
			res = integrate(func = los, a = -1.0*lopmatus, b = 0.0_rk, acc = integreerimise_t2psus, nmin = integration_min_iter, nmax = integration_max_iter)
		else
			res = integrate(func = los, a = 0.0_rk, b = lopmatus, acc = integreerimise_t2psus, nmin = integration_min_iter, nmax = integration_max_iter)
		end if
	contains
		function los(d) result(res)
			implicit none
			real(rk), intent(in) :: d
			real(rk) :: res
			real(rk) :: R, z, theta
			call XcYcl_to_Rztheta(Xc, Yc, d,  comp%sin_incl, comp%cos_incl, comp%tan_incl, comp%sec_incl, comp%theta0, R, z, theta)
			res = ruum_ptr(R,z,theta)
! 			print*, d,res
		end function los
	end function integreerimine_los_dustplane
	
! 	function scalar_over_los(func, default_los_kauguse_piir) result(res)
! 		implicit none
! 		interface
! 			function func(d) result(res)
! 				import rk
! 				real(rk), intent(in) :: d
! 				real(rk) :: res
! 			end function func
! 		end interface
! 		real(rk), intent(in) :: default_los_kauguse_piir
! 		!type(accuracy_type) :: t2psus
! 		real(rk) :: t2psus
! 		real(rk) :: res
! 		real(rk), allocatable, dimension(:,:) :: fvalue
! 		integer :: errcode
! 		t2psus = 1.0e-4*func(0.0_rk)
! ! 		integrate_gauss_adaptive(func,a,b,acc,nmax,nmin,errcode,incracc,fvalue,log_base)
! ! 		res = integrate_gauss_adaptive(func=func, a=-1.0*default_los_kauguse_piir, b=default_los_kauguse_piir, acc=t2psus, nmax=10, nmin=2, errcode=errcode, incracc = .true.)
! 		res = integrate(func, -1.0*default_los_kauguse_piir, default_los_kauguse_piir, t2psus)
! 	end function scalar_over_los
end module los_real_integration_module