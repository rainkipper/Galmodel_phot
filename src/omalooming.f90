module omalooming_module
	!maksimeerib mingit funktsiooni
	!arvatavasti on probleemid kui mingi parameeter l2heneb oma lubatud piiri l2hedale
! 	integer, parameter :: rk = kind(1.0)
	use constants_module
	real(rk), parameter, private :: gold = (3.0-sqrt(5.0))/2.0
	integer, parameter, private :: N_iter_golden = 14
	integer, parameter, private :: N_iter_fittimine = 220
	integer, parameter, private :: mis_versioon = 2
	interface
		function fititav_fun(par_list) result(res)
			import rk
			implicit none
			real(rk), dimension(:), intent(in) :: par_list
			real(rk) :: res
		end function fititav_fun
	end interface
contains
	function fittimine(fun, minpar, maxpar) result(res1)
		implicit none
		procedure(fititav_fun) :: fun
		real(rk), dimension(:), allocatable, intent(in) :: minpar, maxpar !parameetriruumi piirid fittimisel... prior
		real(rk), dimension(:), allocatable :: res1, res2, res3, res4, res5, res6
		real(rk), dimension(1:3) :: valik
		real(rk), dimension(:), allocatable :: X0, X1 !nende punktidega moodustatud sirgel fitib
		real(rk) :: kmin, kmax, kbest !k n2itab kui kaugel X0-st piki sirget X0-X1
		integer :: i, halvim, parim

		allocate(res1(1:size(minpar, 1))); call genereeri_suvaline_X(res1)
		allocate(res2(1:size(minpar, 1))); call genereeri_suvaline_X(res2)
		allocate(res3(1:size(minpar, 1))); call genereeri_suvaline_X(res3)
		allocate(res4(1:size(minpar, 1))); call genereeri_suvaline_X(res4)
		allocate(res5(1:size(minpar, 1))); call genereeri_suvaline_X(res5)
		allocate(res6(1:size(minpar, 1))); call genereeri_suvaline_X(res6)
		allocate(X0(1:size(minpar, 1)))
		allocate(X1(1:size(minpar, 1)))
		call genereeri_suvaline_X(res1) !just random initial value
		
		select case(mis_versioon)
		case(1)
			do i=1,N_iter_fittimine
				print*, i, i*N_iter_golden
				X0 = res1
				call genereeri_suvaline_X(X1)
				call leia_k_piirid(X0, X1, kmin, kmax)
				kmin = kmin * (N_iter_fittimine - i)/(10+N_iter_fittimine) !et kahandab otsitavat ruumi vastavalt progressile
				kmax = kmax * (N_iter_fittimine - i)/(10+N_iter_fittimine) !samaaegselt t2psus suureneb
	! 			print*, kmin, kmax
				kbest = fit_1D(to_1D_fit, kmin, kmax, valik(1))
				res1 = X0 + kbest * (X1 - X0)
	! 			print*, res
			end do
		case(3)
			call genereeri_suvaline_X(res1)
			call genereeri_suvaline_X(res2)
			call genereeri_suvaline_X(res3)
			do i = 1,N_iter_fittimine
				
			end do
		case(2)
			do i=1,N_iter_fittimine
				print*, i, i*N_iter_golden*3
				!1 ja 2 vahel parim
				X0 = res1; X1 = res2
				call leia_k_piirid(X0, X1, kmin, kmax)
				kbest = fit_1D(to_1D_fit, kmin, kmax, valik(1))
				res4 = X0 + kbest * (X1 - X0)

				!1 ja 3 vahel parim
				X0 = res1; X1 = res3
				call leia_k_piirid(X0, X1, kmin, kmax)
				kbest = fit_1D(to_1D_fit, kmin, kmax, valik(2))
				res5 = X0 + kbest * (X1 - X0)
				!2 ja 3 vahel parim
				X0 = res2; X1 = res3
				call leia_k_piirid(X0, X1, kmin, kmax)
				kbest = fit_1D(to_1D_fit, kmin, kmax, valik(3))
				res6 = X0 + kbest * (X1 - X0)

	! 			print*, res
				res1 = res4; res2 = res5; res3 = res6
				halvim = minloc(valik, 1)
				select case(halvim)
					case(1)
						 X0 = 0.5*(res2 + res3);  X1 = res1
					case(2)
						 X0 = 0.5*(res1 + res3);  X1 = res2
					case(3)
						 X0 = 0.5*(res1 + res2);  X1 = res3
				end select
				call leia_k_piirid(X0, X1, kmin, kmax)
				kbest = fit_1D(to_1D_fit, kmin, kmax, valik(halvim))
				select case(halvim)
				case(1); call genereeri_suvaline_X(res1)
				case(2); call genereeri_suvaline_X(res2)
				case(3); call genereeri_suvaline_X(res3)
				end select


			end do
			select case(maxloc(valik, 1))
				case(1); res1 = res4
				case(2); res1 = res5
				case(3); res1 = res6
			end select
	end select
	contains
		function to_1D_fit(k) result(res)
			implicit none
			real(rk), intent(in) :: k
			real(rk) :: res
			real(rk), dimension(:), allocatable :: X
			allocate(X(1:size(X0)))
			X = X0 + k*(X1-X0) !teisendab 1D fittimise parameetrid oige loglike omaks
			res = fun(X)
		end function to_1D_fit
		subroutine genereeri_suvaline_X(res, eelmine)
			implicit none
			real(rk), dimension(:), allocatable, intent(inout) :: res
			real(rk), dimension(:), allocatable, intent(in), optional :: eelmine
			real(rk) :: frac
			real(rk), dimension(1:size(X0)) ::  delta
			integer :: j
			call random_number(res)
			if(present(eelmine)) then
				!variant 1
				!frac = (N_iter_fittimine - i)/real(N_iter_fittimine)
				!res = eelmine*(1.0-frac) + frac*minpar + res*frac*(maxpar - minpar)
				!variant 2
				do j=1,size(X0)
					delta(j) = min( abs(eelmine(j)-minpar(j)), abs(eelmine(j)-maxpar(j)) )
				end do
				res = eelmine + delta  * (res - 0.5)
			else
				res = minpar + res*(maxpar - minpar)
			end if
			
		end subroutine genereeri_suvaline_X
		subroutine leia_k_piirid(X0, X1, kmin, kmax)
			implicit none
			!eesm2rk on leida, milliste k v22rtuste korral on koik parameetrid fittimise piirkonnas sees
			!parameetrite piire kasutab ylemisest funktsioonist
			real(rk), dimension(:), allocatable, intent(in) :: X0, X1
			real(rk), dimension(:), allocatable :: k1, k2
			real(rk), intent(out) :: kmax, kmin
			real(rk) :: tmp
			allocate(k1(1:size(X0))); allocate(k2(1:size(X0))); 
			!leiab koikide parameetrite korral piirid kunas jouab priori servale... 
			k1 = (maxpar-X0)/(X1-X0); k2 = (minpar-X0)/(X1-X0)
			kmax =  minval(k1, 1, k1>0.0); tmp =  minval(k2, 1, k2>0.0); 
			if(tmp<kmax) kmax = tmp !kui teises osas v2iksem, siis vahetab
			kmin =  maxval(k1, 1, k1<0.0); tmp =  maxval(k2, 1, k2<0.0); 
			if(tmp>kmin) kmin = tmp 
		end subroutine leia_k_piirid
	end function fittimine
	
	function fit_1D(fun, min, max, best_funval) result(res)
		!otsib maksimumi 1D funktsioonil
		implicit none
		real(rk), intent(in) :: min, max
		real(rk), intent(out) :: best_funval
		real(rk) :: res
		interface
			function fun(k) result(res)
				import rk
				implicit none
				real(rk), intent(in) :: k
				real(rk) :: res
			end function fun
		end interface
		real(rk), dimension(1:4) :: punktid, f
		integer :: i
		!esimene iteratsioon k2sitsi
		punktid = min + [0.0_rk, gold, 1.0-gold, 1.0_rk]*(max - min)
		do i=1,4
			f(i) = fun(punktid(i))
		end do
		!
		do i=1,N_iter_golden
			if(f(2)>f(3)) then
				!ymber korraldamine ning siis uue v22rtuse arvutamine.
				!selline ymber korraldamine voimaldab koige v2hem uuesti arvutamise punkte.
				punktid(4) = punktid(3)
				punktid(3) = punktid(2)
				punktid(1) = punktid(1)
				punktid(2) = punktid(1)+gold*(punktid(4)-punktid(1))
				f(4) = f(3)
				f(3) = f(2)
				f(1) = f(1)
				f(2) = fun(punktid(2))
			else
				punktid(4) = punktid(4)
				punktid(1) = punktid(2)
				punktid(2) = punktid(3)
				punktid(3) = punktid(1)+(1.0-gold)*(punktid(4)-punktid(1))
				f(4) = f(4)
				f(1) = f(2)
				f(2) = f(3)
				f(3) = fun(punktid(3))
			end if
		end do
		res = punktid(maxloc(f, 1))
		best_funval = maxval(f, 1)
		print*, "1D parim LL = ", maxval(f,1)
	end function fit_1D
end module omalooming_module