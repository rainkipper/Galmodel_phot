module Newton_Rhapson_module
	use constants_module
	use yldine_matemaatika_module
	integer, parameter :: NR_maxiter = 50
contains
	function optim_NR(initial_value, der1, der2) result(res)
		!leiab maksimumi
		implicit none
		real(rk), dimension(:), intent(in) :: initial_value
		interface
			function der1(val) result(res)
				import rk
				implicit none
				real(rk), dimension(:), allocatable, intent(in) :: val
				real(rk), dimension(:), allocatable :: res
			end function der1
			function der2(val) result(res)
				import rk
				implicit none
				real(rk), dimension(:), intent(in) :: val
				real(rk), dimension(:,:), allocatable :: res
			end function der2
		end interface
		real(rk), dimension(:), allocatable :: grad, nihe, res, res_eelmine
		real(rk), dimension(:,:), allocatable :: hessian, inv_hessian
		integer :: iter, N, m, k
		
		N = size(initial_value, 1)
		allocate(res(1:N)); res = initial_value		
		allocate(res_eelmine(1:N)); allocate(nihe(1:N))
		do iter=1,NR_maxiter
			grad = der1(res)
			hessian = der2(res)
			inv_hessian = fun_inv_mx(hessian)
			!
			! ================= edasi liikumine miinimumi poole =================
			!
			inv_hessian = fun_inv_mx(hessian) !poordmaatriksi leidmine, et edasi liikuda
			nihe = 0.0
			res_eelmine = res !salvestab viimase, et testida koonduvust
			do k = 1,N

				do m=1,N
					nihe(k) = nihe(k) + grad(m) * inv_hessian(k,m) !maatriks korrutamine sisuliselt
				end do

				
				nihe(k) = nihe(k) * massi_fittimise_hyppe_kordaja
			end do
			if(any(res - nihe < 0)) then				!kui ekstrapoolib liiga kaugele
				res = res - nihe * minval(abs(massi_fittimise_hyppe_kordaja*res/nihe)) !minval, et l2hima nullile parameetri j2rgi voetaks... siin voib olla erinev massi_fittimise_hyppe_kordaja
			else
				res = res - nihe
			end if
			if(maxval(abs(res-res_eelmine)/res)<massif_fiti_rel_t2psus) then
! 				print*, "t2psus, iter:", minval(abs(res-res_eelmine)/res), iter
				exit
			end if
		end do
	end function optim_NR
end module Newton_Rhapson_module