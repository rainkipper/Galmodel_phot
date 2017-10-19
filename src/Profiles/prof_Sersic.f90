module prof_Sersic_module
	use profiles_den
	type, extends(prof_den_base_type):: prof_Sersic_type
		!sisendparameetrid	
		real(rk) :: M
		real(rk) :: a0
		real(rk) :: N
		real(rk) :: q 
		!lisatakse kiiremaks arvutuseks
		real(rk) :: I0
		real(rk) :: b
		real(rk) :: inv_a0
		real(rk) :: inv_N
		real(rk) :: inv_q
	contains
		procedure :: init_profile 	=> init_prof_Sersic
		procedure :: fun_den 		=> fun_den_Sersic
		procedure :: set_val	     => set_val_Sersic
		procedure :: get_val	    => get_val_Sersic
		procedure :: sanity_check	=> sanity_check_Sersic
		procedure :: fun_los_lopmatus => fun_los_lopmatus_Sersic
	end type prof_Sersic_type
contains
	function fun_los_lopmatus_Sersic(prof, incl, Xc, Yc) result(res)
		class(prof_Sersic_type), intent(in) :: prof
		real(rk), intent(in), optional :: incl, Xc, Yc
		real(rk) :: res
! 		if(present(incl)) then
! 			res = 4*(prof%a0*sin(incl) + prof%a0*prof%q*cos(incl))!
! 		else
			res = default_los_kauguse_piir
! 		end if
	end function fun_los_lopmatus_Sersic
	subroutine set_val_Sersic(prof, par, val)
		implicit none
		class(prof_Sersic_type), intent(inout) :: prof
		character(len=*), intent(in) 	:: par
		real(rk), intent(in) 			:: val
		select case(par)
		case("M"); prof%M = val
		case("a0"); prof%a0 = val; 
		case("N"); prof%N = val
		case("q"); prof%q = val
		case default; stop "Err: niisugust parameetrit bar profiilil pole"
		end select
	end subroutine set_val_Sersic
	subroutine get_val_Sersic(prof, par, val)
		implicit none
		class(prof_Sersic_type), intent(inout) :: prof
		character(len=*), intent(in) 	:: par
		real(rk), intent(out) 			:: val
		select case(par)
		case("M") ;	val = prof%M 
		case("a0");	val = prof%a0
		case("N") ;	val = prof%N 
		case("q") ;	val = prof%q 
		case default; stop "Err: niisugust parameetrit bar profiilil pole"
		end select
		
	end subroutine get_val_Sersic
	elemental function fun_den_Sersic(prof, R,z,theta) result(res)
		implicit none
		class(prof_Sersic_type), intent(in) 	:: prof
		real(rk), intent(in) 				:: R,z
		real(rk), intent(in), optional		:: theta
		real(rk) 							:: res
		real(rk)							:: x,y,a
		x= R*cos(theta)
		y= R*sin(theta)
		a = sqrt(x**2 + (y*prof%inv_q)**2)
		res = prof%I0 * exp(-1.0*prof%b*(a*prof%inv_a0)**prof%inv_N)
	end function fun_den_Sersic
	subroutine init_prof_Sersic(prof)
		implicit none
		class(prof_Sersic_type), intent(inout) :: prof
		!amoeba lolluste vastu... tegelikult vast ei maksaks teha seda:
! 		prof%M = abs(prof%M)
! 		prof%a0 = abs(prof%a0)
! 		prof%N = abs(prof%N)
! 		prof%q = abs(prof%q)
		prof%kas_3D = .false.
		prof%intrinsic_symmetry = "none"
		prof%b = 2*prof%N - 0.324 !voetud Ciotti 1991 toost
		prof%I0 = prof%M/(prof%a0**2*2*pi*prof%N/(prof%b**(2*prof%N)))*gamma(2*prof%N)
		prof%inv_q = 1.0/prof%q
		prof%inv_a0 = 1.0/prof%a0
		prof%inv_N = 1.0/prof%N
	end subroutine init_prof_Sersic
	subroutine sanity_check_Sersic(prof)
		implicit none
		class(prof_Sersic_type), intent(in) :: prof		
		if(prof%M<0) stop("Sanity check for density not passed: bar: M<0")
		if(prof%N<0.5) stop("Sanity check for density not passed: bar: N<0.5... b leidmise prk")
		if(prof%N>10.0) stop("Sanity check for density not passed: bar: N>10... b leidmise prk")
		if(prof%a0<0) stop("Sanity check for density not passed: bar: a0<0")
		if(prof%q<0) stop("Sanity check for density not passed: bar: q<0")
		if(prof%q>1) stop("Sanity check for density not passed: bar: q>1")
	end subroutine sanity_check_Sersic
end module prof_Sersic_module