module prof_Einasto_module
	use profiles_den
	use fgsl

	type, extends(prof_den_base_type):: prof_Einasto_type
		!sisendparameetrid	
		real(rk) :: M
		real(rk) :: a0
		real(rk) :: N
		real(rk) :: q 
		!lisatakse kiiremaks arvutuseks
		real(rk) :: rho0
		real(rk) :: dN
		real(rk) :: k
		real(rk) :: h
		real(rk) :: inv_N
		real(rk) :: inv_ka0
		!tuletise parameetrid
		real(rk) :: der_N1 !abimuutuja
		real(rk) :: psi3N, psi2N
		real(rk) :: N2, N4
		real(rk) :: psi13, psi12
		real(rk) :: q2
! 		procedure(kujuE), pointer :: tuletis
	contains
		procedure :: init_profile 	=> init_prof_Einasto
		procedure :: fun_den 		=> fun_den_Einasto
		procedure :: set_val	     => set_val_Einasto
		procedure :: get_val	    => get_val_Einasto
		procedure :: sanity_check	=> sanity_check_Einasto
		procedure :: fun_los_lopmatus => fun_los_lopmatus_Einasto
		procedure :: lingi_tuletis => lingi_tuletis_Einasto
	end type prof_Einasto_type
! 	procedure(kuju), pointer ::  tuletis_pointer
! 	interface
! 	function kujuE(prof, R,z,theta) result(res)
! 		import rk
! 		import prof_Einasto_type
! 		implicit none
! 		class(prof_Einasto_type), intent(in) 	:: prof
! 		real(rk), intent(in), dimension(:) 				:: R,z
! 		real(rk), intent(in), optional, dimension(:)		:: theta
! 		real(rk), dimension(1:size(R,1)) 							:: res
! 	end function kujuE
! 	end interface

contains
	subroutine lingi_tuletis_Einasto(prof, par1, par2)
		implicit none
		class(prof_Einasto_type), intent(inout) 	:: prof
! 		class(prof_den_base_type), intent(inout) 	:: prof
		character(len=default_character_length), intent(in) :: par1
		character(len=default_character_length), intent(in), optional :: par2
		
		select type(mingiasi=>prof)
		class is (prof_Einasto_type)
			if(present(par2)) then
				if(trim(par1) == "N" .and. trim(par2)=="N") prof%tuletis => den_derivative_NN
				if(trim(par1) == "a0" .and. trim(par2)=="N") prof%tuletis => den_derivative_Na0
				if(trim(par1) == "N" .and. trim(par2)=="a0") prof%tuletis => den_derivative_Na0
				if(trim(par1) == "N" .and. trim(par2)=="q") prof%tuletis => den_derivative_Nq
				if(trim(par1) == "q" .and. trim(par2)=="N") prof%tuletis => den_derivative_Nq
				if(trim(par1) == "a0" .and. trim(par2)=="q") prof%tuletis => den_derivative_a0q
				if(trim(par1) == "q" .and. trim(par2)=="a0") prof%tuletis => den_derivative_a0q
				if(trim(par1) == "q" .and. trim(par2)=="q") prof%tuletis => den_derivative_qq
				if(trim(par1) == "a0" .and. trim(par2)=="a0") prof%tuletis => den_derivative_a0a0
			else
				if(trim(par1) == "N") prof%tuletis => den_derivative_N
				if(trim(par1) == "a0") prof%tuletis => den_derivative_a0
				if(trim(par1) == "q") prof%tuletis => den_derivative_q
			end if
		end select
	end subroutine lingi_tuletis_Einasto
	function den_derivative_N(prof, R, z, theta) result(res)
		implicit none
		class(prof_Einasto_type), intent(in) 	:: prof
		real(rk), intent(in), dimension(:) 				:: R,z
		real(rk), intent(in), optional, dimension(:)		:: theta
		real(rk) 							:: res(1:size(R, 1))
		real(rk)							:: a(1:size(R, 1))
		a = sqrt(R*R + (z/prof%q)**2)
		res = -1.0/(prof%inv_N*prof%inv_N)*(prof%der_N1 - log(a*prof%inv_ka0))*prof%rho0*(a*prof%inv_ka0)**prof%inv_N * exp(-1.0*(a*prof%inv_ka0)**prof%inv_N)
	end function den_derivative_N
	function den_derivative_a0(prof, R, z, theta) result(res)
		implicit none
		class(prof_Einasto_type), intent(in) 	:: prof
		real(rk), intent(in), dimension(:) 				:: R,z
		real(rk), intent(in), optional, dimension(:)		:: theta
		real(rk) 							:: res(1:size(R, 1))
		real(rk)							:: a(1:size(R, 1)), exp_alune(1:size(R, 1)), exp_alune_astmeta(1:size(R, 1))
		a = sqrt(R*R + (z/prof%q)**2)
		exp_alune_astmeta = (a*prof%inv_ka0)
		exp_alune = exp_alune_astmeta**(1/prof%N)
		res = prof%rho0*z**2*exp_alune*exp(-1.0*exp_alune)/((R**2*prof%q**2 + z**2)*prof%N*prof%q)
	end function den_derivative_a0
	function den_derivative_q(prof, R, z, theta) result(res)
		implicit none
		class(prof_Einasto_type), intent(in) 	:: prof
		real(rk), intent(in), dimension(:) 				:: R,z
		real(rk), intent(in), optional, dimension(:)		:: theta
		real(rk) 							:: res(1:size(R, 1))
		real(rk)							:: a(1:size(R, 1)), exp_alune(1:size(R, 1)), exp_alune_astmeta(1:size(R, 1))
		a = sqrt(R*R + (z/prof%q)**2)
		exp_alune_astmeta = (a*prof%inv_ka0)
		exp_alune = exp_alune_astmeta**(1/prof%N)
		res = prof%rho0*exp_alune*exp(-1.0*exp_alune)/(prof%N*prof%a0)
	end function den_derivative_q
	function den_derivative_a0q(prof, R, z, theta) result(res)
		implicit none
		class(prof_Einasto_type), intent(in) 	:: prof
		real(rk), intent(in), dimension(:) 				:: R,z
		real(rk), intent(in), optional, dimension(:)		:: theta
		real(rk) 							:: res(1:size(R, 1))
		real(rk)							:: a(1:size(R, 1)), exp_alune(1:size(R, 1)), exp_alune_astmeta(1:size(R, 1))
		a = sqrt(R*R + (z/prof%q)**2)
		exp_alune_astmeta = (a*prof%inv_ka0)
		exp_alune = exp_alune_astmeta**(1/prof%N)

		!SAGE-ist kopeeritud valem, mida kopeerimise j2rel muditud...
		res = prof%rho0*z**2*(exp_alune - 1)*exp_alune*exp(-1.0*exp_alune)/((R**2*prof%q2 + z*z)*prof%N2*prof%a0*prof%q)
	end function den_derivative_a0q	
	function den_derivative_a0a0(prof, R, z, theta) result(res)
		implicit none
		class(prof_Einasto_type), intent(in) 	:: prof
		real(rk), intent(in), dimension(:) 				:: R,z
		real(rk), intent(in), optional, dimension(:)		:: theta
		real(rk) 							:: res(1:size(R, 1))
		real(rk)							:: a(1:size(R, 1)), exp_alune(1:size(R, 1)), exp_alune_astmeta(1:size(R, 1))
		a = sqrt(R*R + (z/prof%q)**2)
		exp_alune_astmeta = (a*prof%inv_ka0)
		exp_alune = exp_alune_astmeta**(1/prof%N)

		!SAGE-ist kopeeritud valem, mida kopeerimise j2rel muditud...
		res = -1.0*(prof%N - exp_alune + 1)*prof%rho0*exp_alune*exp(-1.0*exp_alune)/(prof%N2*prof%a0**2)
	end function den_derivative_a0a0
	function den_derivative_qq(prof, R, z, theta) result(res)
		implicit none
		class(prof_Einasto_type), intent(in) 	:: prof
		real(rk), intent(in), dimension(:) 				:: R,z
		real(rk), intent(in), optional, dimension(:)		:: theta
		real(rk) 							:: res(1:size(R, 1))
		real(rk)							:: a(1:size(R, 1)), exp_alune(1:size(R, 1)), exp_alune_astmeta(1:size(R, 1))
		a = sqrt(R*R + (z/prof%q)**2)
		exp_alune_astmeta = (a*prof%inv_ka0)
		exp_alune = exp_alune_astmeta**(1/prof%N)

		!SAGE-ist kopeeritud valem, mida kopeerimise j2rel muditud...
		res = -(3*prof%N*R**2*prof%q2 + prof%N*z*z - z*z*(exp_alune + z*z))*prof%rho0*z*z*exp_alune*exp(-1.0*exp_alune)/((R**2*prof%q2 + z*z)**2*prof%N2*prof%q2)
	end function den_derivative_qq
	function den_derivative_Nq(prof, R, z, theta) result(res)
		implicit none
		class(prof_Einasto_type), intent(in) 	:: prof
		real(rk), intent(in), dimension(:) 				:: R,z
		real(rk), intent(in), optional, dimension(:)		:: theta
		real(rk) 							:: res(1:size(R, 1))
		real(rk)							:: a(1:size(R, 1)), exp_alune(1:size(R, 1)), exp_alune_astmeta(1:size(R, 1))
		a = sqrt(R*R + (z/prof%q)**2)
		exp_alune_astmeta = (a*prof%inv_ka0)
		exp_alune = exp_alune_astmeta**(1/prof%N)

		!SAGE-ist kopeeritud valem, mida kopeerimise j2rel muditud...
		res = -1.0*( &
		3*prof%N*exp_alune*prof%psi3N -  &
		2*prof%N*exp_alune*prof%psi2N -  &
		exp_alune*log(exp_alune_astmeta) -  &
		3*prof%N*prof%psi3N + 2*prof%N*prof%psi2N + prof%N +  &
		log(exp_alune_astmeta)  &
		) * prof%rho0*z**2*exp_alune*exp(-1.0*exp_alune)/((R**2*prof%q**2 + z**2)*prof%N**3*prof%q)
	end function den_derivative_Nq
	function den_derivative_Na0(prof, R, z, theta) result(res)
		implicit none
		class(prof_Einasto_type), intent(in) 	:: prof
		real(rk), intent(in), dimension(:) 				:: R,z
		real(rk), intent(in), optional, dimension(:)		:: theta
		real(rk) 							:: res(1:size(R, 1))
		real(rk)							:: a(1:size(R, 1)), exp_alune(1:size(R, 1)), exp_alune_astmeta(1:size(R, 1))
		a = sqrt(R*R + (z/prof%q)**2)
		exp_alune_astmeta = (a*prof%inv_ka0)
		exp_alune = exp_alune_astmeta**(1/prof%N)

		!SAGE-ist kopeeritud valem, mida kopeerimise j2rel muditud...
		res = &
		-1.0*( &
		3*prof%N*exp_alune*prof%psi3N -  &
		2*prof%N*exp_alune*prof%psi2N -  &
		exp_alune*log(exp_alune_astmeta) -  &
		3*prof%N*prof%psi3N + 2*prof%N*prof%psi2N + prof%N +  &
		log(exp_alune_astmeta)) *  &
		prof%rho0*exp_alune*exp(-1.0*exp_alune)/(prof%N**3*prof%a0)
	end function den_derivative_Na0
	function den_derivative_NN(prof, R, z, theta) result(res)
		implicit none
		class(prof_Einasto_type), intent(in) 	:: prof
! 		class(prof_den_base_type), intent(in) 	:: prof
		real(rk), intent(in), dimension(:) 				:: R,z
		real(rk), intent(in), optional, dimension(:)		:: theta
		real(rk) 							:: res(1:size(R, 1))
		real(rk)							:: a(1:size(R, 1)), exp_alune(1:size(R, 1)), exp_alune_astmeta(1:size(R, 1))
		
		a = sqrt(R*R + (z/prof%q)**2)
		exp_alune_astmeta = (a*prof%inv_ka0)
		exp_alune = exp_alune_astmeta**(1/prof%N)

		!SAGE-ist kopeeritud valem, mida kopeerimise j2rel muditud...
		res = ( &
		9*prof%N2*exp_alune*prof%psi3N**2 -  &
		12*prof%N2*exp_alune*prof%psi3N*prof%psi2N +  &
		4*prof%N2*exp_alune*prof%psi2N**2 +  &
		4*prof%N2*prof%N*prof%psi12 -  &
		6*prof%N*exp_alune*log(exp_alune_astmeta)*prof%psi3N -  &
		9*prof%N2*prof%N*prof%psi13 -  &
		9*prof%N2*prof%psi3N**2 +  &
		4*prof%N*exp_alune*log(exp_alune_astmeta)*prof%psi2N +  &
		12*prof%N2*prof%psi3N*prof%psi2N -  &
		4*prof%N2*prof%psi2N**2 + &
		exp_alune*log(exp_alune_astmeta)**2 +  &
		6*prof%N2*prof%psi3N + 6*prof%N*log(exp_alune_astmeta)*prof%psi3N -  &
		4*prof%N2*prof%psi2N - 4*prof%N*log(exp_alune_astmeta)*prof%psi2N -  &
		2*prof%N*log(exp_alune_astmeta) - log(exp_alune_astmeta)**2) * prof%rho0*exp_alune*exp(-exp_alune)/prof%N4
	end function den_derivative_NN
	
	
	
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
		case("a0"); prof%a0 = val; 
		case("N"); prof%N = val
		case("q"); prof%q = val
		case default; stop "Err: niisugust parameetrit Einasto profiilil pole"
		end select
	end subroutine set_val_Einasto
	subroutine get_val_Einasto(prof, par, val)
		implicit none
		class(prof_Einasto_type), intent(inout) :: prof
		character(len=*), intent(in) 	:: par
		real(rk), intent(out) 			:: val
		select case(par)
		case("M") ;	val = prof%M 
		case("a0");	val = prof%a0
		case("N") ;	val = prof%N 
		case("q") ;	val = prof%q 
		case default; stop "Err: niisugust parameetrit Einasto profiilil pole"
		end select
		
	end subroutine get_val_Einasto
	elemental function fun_den_Einasto(prof, R,z,theta) result(res)
		implicit none
		class(prof_Einasto_type), intent(in) 	:: prof
		real(rk), intent(in) 				:: R,z
		real(rk), intent(in), optional		:: theta
		real(rk) 							:: res
		real(rk)							:: a
		if(present(theta) .or. .not.present(theta)) then
			a = sqrt(R**2 + (z/prof%q)**2)
			res =  prof%rho0 * exp(-1.0*(a/(prof%k*prof%a0))**(1/prof%N))
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
! 		prof%dN = 0.0/0.0
! 		prof%M = 1.0_rk
! 		print*, "TODO: edasi t2ita prof"
		prof%inv_ka0 = 1.0/(prof%k * prof%a0)
		prof%der_N1 = 3*prof%N * fgsl_sf_psi(3*prof%N) - 2*prof%N*fgsl_sf_psi(2*prof%N)
		prof%psi3N = fgsl_sf_psi(3*prof%N)
		prof%psi2N = fgsl_sf_psi(2*prof%N)
		prof%N2 = prof%N*prof%N
		prof%N4 = prof%N2*prof%N2
		prof%psi12 = fgsl_sf_psi_n(1, 2*prof%N)
		prof%psi13 = fgsl_sf_psi_n(1, 3*prof%N)
		prof%q2 = prof%q**2
	end subroutine init_prof_Einasto
	subroutine sanity_check_Einasto(prof)
		implicit none
		class(prof_Einasto_type), intent(in) :: prof		
		if(prof%M<0) stop("Sanity check for density not passed: Einasto: M<0")
		if(prof%N<0.5) stop("Sanity check for density not passed: Einasto: N<0.5")
		if(prof%a0<0) stop("Sanity check for density not passed: Einasto: a0<0")
		if(prof%q<0) stop("Sanity check for density not passed: Einasto: q<0")
		if(prof%q>1) stop("Sanity check for density not passed: Einasto: q>1")
	end subroutine sanity_check_Einasto
end module prof_Einasto_module