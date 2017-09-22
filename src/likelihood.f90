module likelihood_module
	use images_module
	use fill_comp_image_module
	use file_operations_module !vaja ainult ajutiselt testimisel
	use konvolutsioon_module
	type masside_arvutamise_tyyp
		!iga vaatluspildi kohta...
		!kasutatakse masside fittimise eristamises muust fittimisest
		real(rk), dimension(:), allocatable :: w !ehk massi komponentide yhikutest ja  ML tulevad kaalud... niipalju kui massi pilte on ka w pikkus
		real(rk), dimension(:,:), allocatable :: I
		real(rk), dimension(:,:), allocatable :: inv_sigma2
		logical, dimension(:,:), allocatable :: mask
		real(rk), dimension(:,:,:), allocatable :: M !sama pikkusega, mis on w... esimen on mis pildi kohta k2ib, ylej22nud on pildi koordinaadid
	end type masside_arvutamise_tyyp
contains
	function calc_log_likelihood(all_comp, images) result(res)
		implicit none
		type(all_comp_type), intent(inout) :: all_comp
		type(image_type), dimension(:), allocatable, intent(in) :: images
		type(comp_image_real_type), dimension(:), allocatable :: mudelid
		real(rk) :: res
		real(rk), dimension(:,:), allocatable :: pilt, pilt_psf
		integer :: i, j, k
		logical :: kas_koik_pildid_samast_vaatlusest !ehk kui koik sama cutout (ruumiliselt ja lahutuselt), siis mudelpilti peab v2he ymber arvutama
		logical :: via_comp_im, kas_los !pildi t2itmise erip2rad
		character(len=default_character_length) :: mida_arvutatakse
		real(rk), dimension(:), allocatable :: weights !ehk M/L suhted fotomeetria korral
		real(rk) :: yhikute_kordaja !10e10Lsun to counts/s
		type(masside_arvutamise_tyyp), dimension(:), allocatable :: to_massfit
		real(rk), dimension(:), allocatable :: lisakaalud_massile
		logical :: kas_fitib_massid_eraldi
		!need peaks tulema mujalt seadetest, mitte k2sitsi
		kas_fitib_massid_eraldi = .true.
		kas_koik_pildid_samast_vaatlusest = .true. 
		via_comp_im = .true.
		kas_los = .true.
		mida_arvutatakse = "Not in use"
		
		!
		! ========== t2psuse leidmine, mida on vaja mudelpildi arvutamiseks=========
		!
		do i=1,all_comp%N_comp
			all_comp%comp(i)%mass_abs_tol = leia_massi_abs_tol(all_comp%comp(i), images)
		end do
		
		!
		! =========== eraldi massifittimise korral asjade initsialiseerimine
		!
		if(kas_fitib_massid_eraldi) then
			allocate(to_massfit(1:size(images,1)))
			do i=1,size(images)
				allocate(to_massfit(i)%w(1:all_comp%N_comp))
				allocate(to_massfit(i)%I(1:size(images(i)%obs,1), size(images(i)%obs,2)))
				allocate(to_massfit(i)%mask(1:size(images(i)%obs,1), size(images(i)%obs,2)))
				allocate(to_massfit(i)%M(1:all_comp%N_comp, 1:size(images(i)%obs,1), size(images(i)%obs,2)))
				allocate(to_massfit(i)%inv_sigma2(1:size(images(i)%obs,1), size(images(i)%obs,2)))
				to_massfit(i)%I = images(i)%obs
				to_massfit(i)%mask = images(i)%mask
				to_massfit(i)%inv_sigma2 = 1.0/(images(i)%sigma**2  + images(i)%sky_noise**2 + abs(images(i)%obs))
			end do
		end if
		
		!
		! ========= komponentide mudelpiltide arvutamised ==============
		!
		if(kas_koik_pildid_samast_vaatlusest) then
			allocate(mudelid(1:all_comp%N_comp))
			do i=1,size(mudelid, 1)
	    		!reaalselt vaja yhe korra ainult teha (koord arvutused sisuslielt)...seega mitteoptimaalsus siin  
				call create_comp_image_from_obs(mudelid(i), images(1))
				call fill_comp_image(all_comp, i, mudelid(i), via_comp_im, kas_los, mida_arvutatakse)
				!kui massid fitib teistest eraldi, siis salvestab massi pildid eraldi
				if(kas_fitib_massid_eraldi) then
					do j=1,size(images)
						to_massfit(j)%M(i,:,:) = mudelid(i)%mx
					end do
				end if
			end do
		else
			print*, "Not yet implemented in calc_log_likelihood"
			stop
		end if
		
		
		!
		! ============= komponentide mudelpiltide kombineerimised ================
		!
		res = 0.0
		allocate(weights(1:all_comp%N_comp))
		do i=1,size(images)			
			!
			! ================== kaalude arvutamised ========= 
			!
			weights = -1.234 !ehk algselt koigile mingi tobe v22rtus, et hiljem saaks vigasust kontrollida
			do j=1,size(weights, 1)
					weights(j) = images(i)%filter%mass_to_obs_unit(all_comp%comp(1)%dist, all_comp%comp(j)%population_name)
			end do
			if(any(weights == -1.234)) then
				print*, "Vale M/L sisend"; stop
			end if
			if(kas_fitib_massid_eraldi) to_massfit(i)%w(:) = weights(:) !hilisemaks kasutamiseks kui otsustab eraldi fittida
		end do
			
		
		!
		! ========= mudelpiltide kokkupanek ===========
		!		
		if(kas_fitib_massid_eraldi) then
			!lisab kaaludele veel massi kordaja sisemisest fittimisest
			lisakaalud_massile = fiti_massi_kordajad(to_massfit)
			do i=1,size(weights)
				weights(i) = weights(i) * lisakaalud_massile(i)
			end do
		end if
		
		do i=1,size(images)	
			if(allocated(pilt)) deallocate(pilt)
			if(kas_koik_pildid_samast_vaatlusest) then
				pilt = combine_comp_images_to_make_image(mudelid, weights)
			else
				print*, "Not yet implemented in likelihood"
				stop
			end if
			
			!
			! ========== psf rakendamine =========
			!
			if(.false.) then
				if(allocated(pilt_psf)) deallocate(pilt_psf)
				call convolve(pilt, images(i)%psf, pilt_psf)
				print*, "konvoleeritud"
			else
				if(allocated(pilt_psf)) deallocate(pilt_psf)
				allocate(pilt_psf(1:size(pilt,1), 1:size(pilt,2)))
				pilt_psf = pilt
			end if
! 			print*, size(pilt_psf, 1), size(pilt_psf, 2)
			call write_matrix_to_fits(pilt_psf, images(i)%output_mdl_file)
			
			!
			! ======== loglike ise ========
			!
			res = res + sum(-1.0*( (pilt_psf-images(i)%obs)**2*0.5/((images(i)%sigma)**2  + (images(i)%sky_noise**2 + abs(images(i)%obs)))), images(i)%mask) 
		end do
		
		print*, "LL = ", res
	end function calc_log_likelihood
	function leia_massi_abs_tol(comp, images) result(res)
		!leiab kauguste, filtrite jm pohjal massi hajuvuse ning korrutab konstandiga, et saada hajumisest t2psus
		implicit none
		type(comp_type), intent(in) :: comp
		type(image_type), intent(in), dimension(:), allocatable :: images
		real(rk) :: res
		real(rk) :: tmp
		integer :: i
		res = 1.0e10 !suvaline suur suurus algseks
		do i=1,size(images, 1)
			tmp = 1.0/images(i)%filter%mass_to_obs_unit(comp%dist, comp%population_name) !poordvaartus, kuna vaja massi hajumist saada countside hajumisest
			if(tmp<res) then
				res = tmp
			end if
		end do
		res = res * massi_abs_tol_kordaja
		
	end function leia_massi_abs_tol
	function fiti_massi_kordajad_lin_regressioon(to_massfit) result(res)
		!tekitab negatiivseid massi kordajaid, ehk kasutada ainult väga mittekõdunud profiilide korral
		implicit none
		type(masside_arvutamise_tyyp), dimension(:), allocatable :: to_massfit
		integer :: i,j,k
		integer :: N_k
		real(rk), dimension(:), allocatable :: res
		real(rk), dimension(:,:), allocatable :: A, inv_A
		real(rk), dimension(:), allocatable :: B
		stop "Liiga kontrollimata, et toimiks"
		N_k = size(to_massfit(1)%w, 1)
		allocate(B(1:N_k))
		allocate(A(1:N_k, 1:N_k)); allocate(inv_A(1:N_k, 1:N_k))
		do k=1,N_k
			!B leidmine
			B(k) = 0
			do i=1,size(to_massfit)
				B(k) = B(k) + sum(  (to_massfit(i)%I*to_massfit(i)%w(k)*to_massfit(i)%M(k,:,:))*to_massfit(i)%inv_sigma2 , to_massfit(i)%mask )
			end do
			!A leidmine
			do j = 1,N_k
! 			A(k,j)  = 0.0
			A(j,k)  = 0.0
				do i=1,size(to_massfit)
! 					A(k,j)  = A(k,j) + &
					A(j,k)  = A(j,k) + &
					sum( to_massfit(i)%w(k) * to_massfit(i)%M(k,:,:) * to_massfit(i)%w(j) * to_massfit(i)%M(j,:,:) * to_massfit(i)%inv_sigma2, to_massfit(i)%mask )
				end do
			end do
		end do
		!poordmaatriksi leidmine ja massile kaalude saamine
		allocate(res(1:N_k))
		inv_A = fun_inv_mx(A)
		do k=1,N_k
			res(k) = sum( inv_A(k,:) * B(:) )
		end do
		
		print "(5F)", res
		res = (res)
	end function fiti_massi_kordajad_lin_regressioon
	function fiti_massi_kordajad(to_massfit) result(res)
		implicit none
		type(masside_arvutamise_tyyp), dimension(:), allocatable :: to_massfit
		integer :: i,j,k,m, iter
		integer :: N_k, N_i, N_iter
		real(rk), dimension(:), allocatable :: res
		real(rk), dimension(:,:), allocatable :: L_km, L0_km, B_km, inv_L_km !kogu chisq, piltidest chisq, barrier osa chisq-st, poordmaatriks
		real(rk), dimension(:), allocatable :: L_k, L0_k, B_k, nihe
		logical :: kas_barrier
		real(rk), dimension(:,:), allocatable :: tmp_pilt
		real(rk) :: lambda, gamma
real(rk) :: t1,t2, tt1, tt2, dt_grad, dt_hess, dt_nihe
integer :: vidin
		kas_barrier = .true.
		lambda = 50.0 !m22rab kui t2pselt ei tohi massid nulli minna... voib olla problemaatiline kui on suured hypped iteratsioonide vahel.. 
		gamma = 0.7 !ehk kui kiiresti liigub iteratsioonide vahel
		N_iter = 15
		!initsialiseerimised
		N_k = size(to_massfit(1)%w, 1) !komponentide arv
		N_i = size(to_massfit, 1) !piltide arv
		allocate(L0_k(1:N_k)) 
		allocate(L_k(1:N_k)) 
		allocate(L0_km(1:N_k, 1:N_k)); 
		allocate(L_km(1:N_k, 1:N_k))
		if(kas_barrier) then
			allocate(B_km(1:N_k, 1:N_k)); 
			allocate(B_k(1:N_k))
		end if
		allocate(nihe(1:N_k))
		allocate(res(1:N_k)); res = 1.0 !algne masside kordajad on 1.0... seal hakkab edasi roomama
! call  cpu_time(t1); dt_grad=0.0; dt_hess = 0.0
! do vidin = 1,100
! do i=1,N_k; call random_number(res(i)) ;end do; res = res+0.1
! print "(A,5F10.5)", "Enne",res
		do iter = 1, N_iter
			!iteratsiooni ettevalmistus
			L0_k = 0.0; L_k = 0.0; L0_km = 0.0 ; L_km = 0.0 
			if(kas_barrier) then
				B_km = 0.0;  B_k = 0.0
			end if
			nihe = 0.0
			!Hessiani ja gradiendi arvuamised
			do k = 1,N_k
				!
				! ================= gradiendi arvutamine ================= 
				!
! call  cpu_time(tt1)
				do i=1,N_i
					if(allocated(tmp_pilt)) deallocate(tmp_pilt)
					allocate(tmp_pilt(1:size(to_massfit(i)%I, 1), 1:size(to_massfit(i)%I, 2))); tmp_pilt = 0.0
					!selle iteratsiooni heleduse pilt
					do j=1,N_k
						tmp_pilt = tmp_pilt + res(j)*to_massfit(i)%w(j)*to_massfit(i)%M(j,:,:)
					end do
					!likelihoodi gradienti lisamine 
					L0_k(k) = L0_k(k) + sum(2.0*to_massfit(i)%inv_sigma2 * to_massfit(i)%w(k)*to_massfit(i)%M(k,:,:)*(tmp_pilt - to_massfit(i)%I), to_massfit(i)%mask)
				end do
				L_k(k) = L0_k(k)
				if(kas_barrier) L_k(k) = L_k(k) - lambda / res(k) !kui piirab masse seestpoolt
! call  cpu_time(tt2)
! dt_grad = dt_grad + tt2-tt1
				!
				! ================= Hessiani komopnendid ================= 
				!
				do m=k,N_k
					do i=1,N_i
						L0_km(k,m) = L0_km(k,m) + sum(2*to_massfit(i)%inv_sigma2 * to_massfit(i)%w(k)*to_massfit(i)%M(k,:,:)*to_massfit(i)%w(m)*to_massfit(i)%M(m,:,:) ,to_massfit(i)%mask)
					end do
					L_km(k,m) = L0_km(k,m)
					L_km(m,k) = L_km(k,m) !symmeetrilise maatriksi t2itmine
				end do
				if(kas_barrier) L_km(k,k) = L_km(k,k) + lambda/res(k)**2 !kui piirab masse seestpoolt... sisaldab ainult diagonaalil olevaid elemente
! call  cpu_time(tt1)
! dt_hess = dt_hess - tt2+tt1
			end do

			!
			! ================= edasi liikumine miinimumi poole =================
			!
! call  cpu_time(tt1)
			inv_L_km = fun_inv_mx(L_km) !poordmaatriksi leidmine, et edasi liikuda
			nihe = 0.0
			do k = 1,N_k
				do m=1,N_k
					nihe(k) = nihe(k) + L_k(m) * inv_L_km(k,m) !maatriks korrutamine sisuliselt
				end do
				nihe(k) = nihe(k) * gamma
			end do
			if(any(res - nihe < 0)) then
				!kui ekstrapoolib liiga kaugele
				res = res - nihe * minval(abs(0.9*res/nihe)) !minval, et l2hima nullile parameetri j2rgi voetaks
			else
				res = res - nihe
			end if
! print "(5F10.5)", res
! call  cpu_time(tt2)
! dt_nihe = tt2-tt1
		end do
! print "(5F10.5)", res
! end do
! call  cpu_time(t2)
! print*, "fitting done:D", t2-t1
! print "(A,3F10.5)", "fitting done:D", dt_grad, dt_hess, dt_nihe
	end function fiti_massi_kordajad
	
end module likelihood_module