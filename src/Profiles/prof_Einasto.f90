module prof_Einasto_module
	use profiles_den
	type, extends(prof_den_base_type):: prof_Einasto_type
		!sisendparameetrid	
		real(rk) :: M, logM
		real(rk) :: a0
		real(rk) :: N
		real(rk) :: q 
		!lisatakse kiiremaks arvutuseks
		real(rk) :: rho0
		real(rk) :: dN
		real(rk) :: k
		real(rk) :: h
		real(rk) :: inv_N
		real(rk) :: inv_q2
		real(rk) :: inv_2N
		real(rk) :: x0
	contains
		procedure :: init_profile 	=> init_prof_Einasto
		procedure :: fun_den 		=> fun_den_Einasto
		procedure :: set_val	     => set_val_Einasto
		procedure :: get_val	    => get_val_Einasto
		procedure :: sanity_check	=> sanity_check_Einasto
		procedure :: fun_los_lopmatus => fun_los_lopmatus_Einasto
	end type prof_Einasto_type
contains
	function fun_los_lopmatus_Einasto(prof, incl, Xc, Yc) result(res)
		class(prof_Einasto_type), intent(in) :: prof
		real(rk), intent(in), optional :: incl, Xc, Yc
		real(rk) :: res
		if(present(incl)) then
			res = 4*(prof%a0*sin(incl) + prof%a0*prof%q*cos(incl))!
		else
			res = default_los_kauguse_piir
		end if
	end function fun_los_lopmatus_Einasto
	subroutine set_val_Einasto(prof, par, val)
		implicit none
		class(prof_Einasto_type), intent(inout) :: prof
		character(len=*), intent(in) 	:: par
		real(rk), intent(in) 			:: val
		select case(par)
		case("M"); prof%M = val
		case("logM"); prof%logM = val; prof%M = 10**prof%logM
		case("a0"); prof%a0 = val; 
		case("N"); prof%N = val
		case("q"); prof%q = val
		case default; 
			print*, "Err: niisugust parameetrit Einasto profiilil pole:", trim(par)
			stop
		end select
	end subroutine set_val_Einasto
	subroutine get_val_Einasto(prof, par, val)
		implicit none
		class(prof_Einasto_type), intent(in) :: prof
		character(len=*), intent(in) 	:: par
		real(rk), intent(out) 			:: val
		select case(par)
		case("M") ;	val = prof%M 
		case("logM") ;	val = prof%logM 
		case("a0");	val = prof%a0
		case("N") ;	val = prof%N 
		case("q") ;	val = prof%q 
		case default; 
		print*, "Err: niisugust parameetrit Einasto profiilil pole:", trim(par)
		stop
		end select
		
	end subroutine get_val_Einasto
	elemental function fun_den_Einasto(prof, R,z,theta) result(res)
		implicit none
		class(prof_Einasto_type), intent(in) 	:: prof
		real(rk), intent(in) 				:: R,z
		real(rk), intent(in), optional		:: theta
		real(rk) 							:: res
		real(rk)							:: a, x
		if(present(theta) .or. .not.present(theta)) then
! 			a = sqrt(R**2 + (z/prof%q)**2)
! 			res =  prof%rho0 * exp(-1.0*(-a*prof%inv_ka0)**prof%inv_2N)
			x = (R*R + z*z*prof%inv_q2)**prof%inv_2N
			res = prof%rho0 * exp(x * prof%x0)
 
		end if
	end function fun_den_Einasto
	subroutine init_prof_Einasto(prof)
		implicit none
		class(prof_Einasto_type), intent(inout) :: prof
		
		prof%kas_3D = .true.
		!amoeba lolluste vastu... tegelikult vast ei maksaks teha seda:
		prof%M = abs(prof%M)
		prof%a0 = abs(prof%a0)
		prof%N = abs(prof%N)
		prof%q = abs(prof%q)
		
		
		!
		prof%intrinsic_symmetry = "Einasto"
		prof%h = gamma(3*prof%N)**2 / (prof%N * gamma(2*prof%N)**3)
		prof%k = gamma(2*prof%N)/gamma(3*prof%N)
		prof%rho0 = prof%h * prof%M / (4*pi*prof%q * prof%a0**3)
		prof%inv_q2 = 1.0/(prof%q**2)
		prof%inv_2N = 0.5/prof%N
		prof%x0 = -1.0*(1.0/(prof%k * prof%a0))**(1.0/prof%N)
! 		prof%dN = 0.0/0.0
! 		prof%M = 1.0_rk
! 		print*, "TODO: edasi t2ita prof"
	end subroutine init_prof_Einasto
	subroutine sanity_check_Einasto(prof)
		implicit none
		class(prof_Einasto_type), intent(in) :: prof		
		if(prof%M<0) stop("Sanity check for density not passed: Einasto: M<0")
		if(prof%N<0.2) stop("Sanity check for density not passed: Einasto: N<0.5")
		if(prof%a0<0) stop("Sanity check for density not passed: Einasto: a0<0")
		if(prof%q<0) stop("Sanity check for density not passed: Einasto: q<0")
		if(prof%q>1) stop("Sanity check for density not passed: Einasto: q>1")
	end subroutine sanity_check_Einasto
end module prof_Einasto_module