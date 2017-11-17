module prof_dustplane_module
	use profiles_den
	type, extends(prof_den_base_type):: prof_dustplane_type
		!sisendparameetrid	
		real(rk) :: M
		real(rk) :: tau0
	contains
		procedure :: init_profile 	=> init_prof_dustplane
		procedure :: fun_den 		=> fun_den_dustplane
		procedure :: set_val	     => set_val_dustplane
		procedure :: get_val	    => get_val_dustplane
		procedure :: sanity_check	=> sanity_check_dustplane
		procedure :: fun_los_lopmatus => fun_los_lopmatus_dustplane
	end type prof_dustplane_type
contains
	function fun_los_lopmatus_dustplane(prof, incl, Xc, Yc) result(res)
		class(prof_dustplane_type), intent(in) :: prof
		real(rk), intent(in), optional :: incl, Xc, Yc
		real(rk) :: res
		res = default_los_kauguse_piir
	end function fun_los_lopmatus_dustplane
	subroutine set_val_dustplane(prof, par, val)
		implicit none
		class(prof_dustplane_type), intent(inout) :: prof
		character(len=*), intent(in) 	:: par
		real(rk), intent(in) 			:: val
		select case(par)
		case("M"); prof%M = val
		case("tau0"); prof%tau0 = val; 
		case default; stop "Err: niisugust parameetrit dustplane profiilil pole"
		end select
	end subroutine set_val_dustplane
	subroutine get_val_dustplane(prof, par, val)
		implicit none
		class(prof_dustplane_type), intent(in) :: prof
		character(len=*), intent(in) 	:: par
		real(rk), intent(out) 			:: val
		select case(par)
		case("M") ;	val = prof%M 
		case("tau0");	val = prof%tau0
		case default; stop "Err: niisugust parameetrit dustplane profiilil pole"
		end select
		
	end subroutine get_val_dustplane
	elemental function fun_den_dustplane(prof, R,z,theta) result(res)
		implicit none
		class(prof_dustplane_type), intent(in) 	:: prof
		real(rk), intent(in) 				:: R,z
		real(rk), intent(in), optional		:: theta
		real(rk) 							:: res
		real(rk)							:: a
		res = (R-R)/(R-R)
	end function fun_den_dustplane
	subroutine init_prof_dustplane(prof)
		implicit none
		class(prof_dustplane_type), intent(inout) :: prof
		
		prof%kas_3D = .true.
		!amoeba lolluste vastu... tegelikult vast ei maksaks teha seda:
		prof%M = abs(prof%M)
	end subroutine init_prof_dustplane
	subroutine sanity_check_dustplane(prof)
		implicit none
		class(prof_dustplane_type), intent(in) :: prof		
	end subroutine sanity_check_dustplane
end module prof_dustplane_module