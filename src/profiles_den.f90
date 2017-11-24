module profiles_den
	use constants_module
	 
	
	type :: prof_den_base_type !this will be extended to more suitable profile ... e.g. Einasto
		logical :: kas_3D = .false. !ehk kas tuleb yle vaatejoone integreerida
		character(len=default_character_length) :: den_prof_name !
		character(len=default_character_length) :: intrinsic_symmetry = "none"  
	contains
		procedure :: init_profile 	=> init_prof_default
		procedure :: fun_den 		=> fun_den_default
		procedure :: set_val	    => set_val_default
		procedure :: get_val	    => get_val_default
		procedure :: sanity_check	=> sanity_check_default
		procedure :: fun_los_lopmatus => fun_los_lopmatus_default
	end type prof_den_base_type
	
	


contains
	function fun_default(prof, R,z,theta) result(res)
		class(prof_den_base_type), intent(in) :: prof
		real(rk), intent(in) :: R
		real(rk), intent(in) :: z
		real(rk), intent(in), optional :: theta
		real(rk) :: res
		res = 0.0/0.0
		stop "Err: default does not suit for any profile"
		print*, len(prof%den_prof_name), R,z,theta
	end function fun_default
	function fun_los_lopmatus_default(prof, incl, Xc, Yc) result(res)
		class(prof_den_base_type), intent(in) :: prof
		real(rk), intent(in), optional :: incl, Xc, Yc
		real(rk) :: res
		res = default_los_kauguse_piir
	end function fun_los_lopmatus_default
	subroutine set_val_default(prof, par, val)
		implicit none
		class(prof_den_base_type), intent(inout) :: prof
		character(len=*), intent(in) 	:: par
		real(rk), intent(in) 			:: val
		stop "Err: default does not suit for any profile"
		print*, val, len(prof%den_prof_name), par
	end subroutine set_val_default
	subroutine get_val_default(prof, par, val)
		implicit none
		class(prof_den_base_type), intent(in) :: prof
		character(len=*), intent(in) 	:: par
		real(rk), intent(out) 			:: val
		val = 0.0_rk
		stop "Err: default does not suit for any profile"
		print*, len(prof%den_prof_name), par
	end subroutine get_val_default
	elemental function fun_den_default(prof, R,z,theta) result(res)
		implicit none
		class(prof_den_base_type), intent(in) 		:: prof
		real(rk), intent(in) 						:: R,z
		real(rk), intent(in), optional				:: theta
		real(rk)  									:: res
		
		res = 0/0 + R+z+theta + len(prof%den_prof_name)
		
	end function fun_den_default
	subroutine init_prof_default(prof)
		implicit none
		class(prof_den_base_type), intent(inout) :: prof
		stop "Err: default does not suit for any profile"
		print*, len(prof%den_prof_name)
	end subroutine init_prof_default
	subroutine sanity_check_default(prof)
		implicit none
		class(prof_den_base_type), intent(in) :: prof
		logical :: kas_ok
		kas_ok = .false.
		if(kas_ok) stop("Sanity check not passed by default den profile")
		print*, len(prof%den_prof_name)
	end subroutine sanity_check_default
	
	

	
end module profiles_den