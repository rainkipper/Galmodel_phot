module approx_Einasto_module
	use constants_module
	type :: approx_Einasto_type
		real(rk) :: logx0, logx1, logx2
		real(rk) :: f0,f1,f2
		real(rk) :: Q_obs2, q !inverse n2iv**2 ja tegelik lapikus... t2idetakse ainult pealmisel kihil
		integer :: level
		logical :: last_level
		type(approx_Einasto_type), pointer :: next_left, next_right
		end type approx_Einasto_type
		type(approx_Einasto_type), private :: prof !privaatne levinud nime tottu
! 		integer, save :: mitu_E
contains
	subroutine fill_approx_Einasto(los_val_func, incl, q)
		implicit none
		real(rk), intent(in) :: incl, q !n2iv ja tegelik lapikus
		interface
			function los_val_func(Xc, Yc) result(res)
				import rk
				real(rk) :: res
				real(rk), intent(in) :: Xc, Yc
			end function los_val_func
		end interface
! 		mitu_E = 0
		prof%q = q
		prof%Q_obs2 = cos(incl)*cos(incl); prof%Q_obs2 = 1.0/(prof%Q_obs2 + (1.0-prof%Q_obs2)*q*q)
		prof%logx0 = (0.001_rk)
		prof%logx2 = (200.0_rk)
! 		prof%logx0 = log10(0.001_rk)
! 		prof%logx2 = log10(adaptive_image_x1_default)
		prof%logx1 = 0.5*(prof%logx0 + prof%logx2)
		prof%f0 = los_val_func(prof%logx0, 0.0_rk)
		prof%f1 = los_val_func(prof%logx1, 0.0_rk)
		prof%f2 = los_val_func(prof%logx2, 0.0_rk)
! 		prof%f0 = los_val_func(10**prof%logx0, 0.0_rk)
! 		prof%f1 = los_val_func(10**prof%logx1, 0.0_rk)
! 		prof%f2 = los_val_func(10**prof%logx2, 0.0_rk)
		
		prof%last_level = .false.
		prof%level = 0
		call fill_single(prof, prof%logx0, prof%logx2, prof%f0, prof%f2, prof%level+1)
! 		print*, "E arvutusi", mitu_E
	contains
		recursive subroutine fill_single(cell, logx0, logx2, f0, f2, level)
			implicit none
			type(approx_Einasto_type), intent(inout) :: cell
			real(rk) :: l2hend
			real(rk), intent(in) :: logx0, logx2, f0, f2
			integer, intent(in) :: level
! 			mitu_E = mitu_E + 1
			cell%level = level
			cell%logx0 = logx0
			cell%logx2 = logx2
			cell%f0 = f0
			cell%f2 = f2
			cell%logx1 = 0.5*(cell%logx0 + cell%logx2)
! 			cell%f1 = los_val_func(10**cell%logx1, 0.0_rk)
			cell%f1 = los_val_func(cell%logx1, 0.0_rk)
! 			l2hend = (cell%f2 - cell%f0)/(cell%logx2-cell%logx0) * (cell%logx1 - cell%logx0) + cell%f0
			l2hend = 0.5*(cell%f2 + cell%f0)
			!test kas piisavalt arvutatud
			cell%last_level = abs(l2hend - cell%f1)/min(abs(l2hend), abs(cell%f1)) < adaptive_image_edasijagamise_threshold
			if(level < adaptive_image_minlevel) cell%last_level = .false.
			if(level > adaptive_image_maxlevel) cell%last_level = .true.
			if(cell%last_level) then
				return
			else
				allocate(cell%next_right)
				call fill_single(cell%next_right, cell%logx1, cell%logx2, cell%f1, cell%f2, level + 1)
				allocate(cell%next_left); 
				call fill_single(cell%next_left, cell%logx0, cell%logx1, cell%f0, cell%f1, level + 1)
			end if
		end subroutine fill_single
	end subroutine fill_approx_Einasto
	subroutine delete_approx_Einasto()
		implicit none
		call delete_single(prof%next_right)
		deallocate(prof%next_right)
		call delete_single(prof%next_left)
		deallocate(prof%next_left)
	contains
		recursive subroutine delete_single(cell)
			implicit none
			type(approx_Einasto_type), intent(inout) :: cell
			if(.not.cell%last_level) then
				call delete_single(cell%next_left)
				call delete_single(cell%next_right)
				deallocate(cell%next_left)
				deallocate(cell%next_right)
			end if
		end subroutine delete_single
	end subroutine delete_approx_Einasto
	function otsi_approx_Einasto(Xc, Yc) result(res)
		implicit none
		real(rk), intent(in) :: Xc, Yc
		real(rk) :: logA
		real(rk) :: res
		logA = sqrt(Xc*Xc + Yc*Yc*prof%Q_obs2) 
		logA = (logA)
! 		logA = log10(logA)
		if(logA < prof%logx0 .or. logA>prof%logx2) then
! 			print*, "Einasto prof approx piiridest v2ljas", 10**logA
			res = 0.0_rk
		else
			res = otsi_yks(prof, logA)
		end if
		
	contains
		recursive function otsi_yks(cell, logA) result(res)
			implicit none
			type(approx_Einasto_type), intent(in) :: cell
			real(rk), intent(in) :: logA
			real(rk) :: res
			if(.not. cell%last_level) then
				if(logA < cell%logx1) then
					res = otsi_yks(cell%next_left, logA)
				else
					res = otsi_yks(cell%next_right, logA)
				end if
			else
				!v1 = Lagrange polynomial
	! 				res = 0.0
	! 				res = res + (logA - )
				!v2 = linear interpolation kahest osast
				if(logA<cell%logx1) then
					res = (cell%f1 - cell%f0)/(cell%logx1 - cell%logx0)*(logA - cell%logx0) + cell%f0
				else
					res = (cell%f2 - cell%f1)/(cell%logx2 - cell%logx1)*(logA - cell%logx1) + cell%f1
				end if
				!v3  = lineaarne j2medal vorel
! 				res = (cell%f2 - cell%f0)/(cell%logx2 - cell%logx0) * (logA-cell%logx0) + cell%f0
				
			end if
		end function otsi_yks
	end function otsi_approx_Einasto
end module approx_Einasto_module