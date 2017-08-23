module profiles_den
	use constants_module
	
	
	type :: prof_den_base_type !this will be extended to more suitable profile ... e.g. Einasto
		character(len=default_character_length) :: den_prof_name !
		character(len=default_character_length) :: intrinsic_symmetry = "none"  
	contains
		procedure :: init_profile 	=> init_prof_default
		procedure :: fun_den 		=> fun_den_default
		procedure :: set_val	    => set_val_default
		procedure :: get_val	    => get_val_default
		procedure :: sanity_check	=> sanity_check_default
	end type prof_den_base_type
	
	


contains
	function fun_default(prof, R,z,theta) result(res)
		class(prof_den_base_type), intent(in) :: prof
		real(rk), intent(in) :: R
		real(rk), intent(in) :: z
		real(rk), intent(in), optional :: theta
		real(rk) :: res
		res = 0.0/0.0
	end function fun_default
	subroutine set_val_default(prof, par, val)
		implicit none
		class(prof_den_base_type), intent(inout) :: prof
		character(len=*), intent(in) 	:: par
		real(rk), intent(in) 			:: val
		stop "Err: default does not suit for any profile"
	end subroutine set_val_default
	subroutine get_val_default(prof, par, val)
		implicit none
		class(prof_den_base_type), intent(inout) :: prof
		character(len=*), intent(in) 	:: par
		real(rk), intent(out) 			:: val
		val = 0.0_rk
		stop "Err: default does not suit for any profile"
	end subroutine get_val_default
	elemental function fun_den_default(prof, R,z,theta) result(res)
		implicit none
		class(prof_den_base_type), intent(in) 		:: prof
		real(rk), intent(in) 						:: R,z
		real(rk), intent(in), optional				:: theta
		real(rk)  									:: res
		res = 0/0
	end function fun_den_default
	subroutine init_prof_default(prof)
		implicit none
		class(prof_den_base_type), intent(inout) :: prof
		stop "Err: default does not suit for any profile"
	end subroutine init_prof_default
	subroutine sanity_check_default(prof)
		implicit none
		class(prof_den_base_type), intent(in) :: prof
		logical :: kas_ok
		kas_ok = .false.
		if(kas_ok) stop("Sanity check not passed by default den profile")
	end subroutine sanity_check_default
	
	

	
end module profiles_den