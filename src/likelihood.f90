module likelihood_module
	use images_module
	use fill_comp_image_module
	use file_operations_module !vaja ainult ajutiselt testimisel
	use psf_rakendamine_module
	type masside_arvutamise_tyyp
		!iga vaatluspildi kohta...
		!kasutatakse masside fittimise eristamises muust fittimisest
		real(rk), dimension(:), allocatable :: w !ehk massi komponentide yhikutest ja  ML tulevad kaalud... niipalju kui massi pilte on ka w pikkus
		real(rk), dimension(:,:), allocatable :: I !pilt ise
		real(rk), dimension(:,:), allocatable :: inv_sigma2 !kaalud
		logical, dimension(:,:), allocatable :: mask !
		real(rk), dimension(:,:,:), allocatable :: M !sama pikkusega, mis on w... esimen on mis pildi kohta k2ib, ylej22nud on pildi koordinaadid
	end type masside_arvutamise_tyyp
	!globaalsed muutujad selle mooduli jaoks
	type(comp_image_real_type), dimension(:), allocatable, private :: mudelid  !siin hoitakse mudeli andmeid ... init_log_likelihood juures pannakse paika
	type(masside_arvutamise_tyyp), dimension(:), allocatable, private :: to_massfit !lihtsustav muutuja
	logical, parameter, private :: kas_fitib_massid_eraldi = .true.
	logical, parameter, private :: kas_koik_pildid_samast_vaatlusest = .true. 
	logical, parameter, private :: via_adaptive_im = .false.
	logical, parameter, private :: kas_los = .false.
	logical, parameter, private :: kas_barrier = .true.
	logical, parameter, private :: kas_rakendab_psf = .true.
	real(rk), parameter,private :: massif_fiti_rel_t2psus = 0.003 !suhteline t2psus, mille korral loeb koondunuks masside eraldi fittimise
	real(rk), dimension(:), allocatable, private :: massi_kordajad_eelmine !massi kordajad... globaalne muutuja, et j2rgmine loglike arvutamine oleks hea algl2hend votta
	integer, save :: LL_counter = 0 !lihtsalt, mitu LL juba arvutatud
contains
	subroutine init_calc_log_likelihood(all_comp, images)
		!eesm2rk on teha m2llu ja initsialiseerida piltide arvuamise asjad
		implicit none
		type(all_comp_type), intent(inout) :: all_comp
		type(image_type), dimension(:), allocatable, intent(in) :: images
		integer :: i
		!
		! =========== eraldi massifittimise korral asjade initsialiseerimine
		!
		allocate(mudelid(1:all_comp%N_comp))
		if(kas_koik_pildid_samast_vaatlusest) then
			do i=1,size(mudelid)
				call create_comp_image_from_obs(mudelid(i), images(1))
			end do
		else
			stop "not ready in likelihood"
		end if
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
		else
			stop "not implemented in init log likelihood"
		end if

		
		mudelid(:)%recalc_image = .true.
	end subroutine init_calc_log_likelihood
	function calc_log_likelihood(all_comp, images, lisakaalud_massile) result(res)
		implicit none
		type(all_comp_type), intent(inout) :: all_comp
		type(image_type), dimension(:), allocatable, intent(in) :: images
		real(rk) :: res
		real(rk), dimension(:,:), allocatable :: pilt, pilt_psf
		integer :: i, j !, k
		character(len=default_character_length) :: mida_arvutatakse
		real(rk), dimension(:,:), allocatable :: weights !ehk M/L suhted fotomeetria korral ... esimene indeks pilt, teine komponent
		real(rk), dimension(:), allocatable :: weights_for_single_im
		real(rk), dimension(:), allocatable, intent(out) :: lisakaalud_massile !optional output
		!need peaks tulema mujalt seadetest, mitte k2sitsi

call testi_psf(images(1)%obs, images(1)%psf)
		mida_arvutatakse = "Not in use"
		LL_counter = LL_counter + 1
		!
		! ========== t2psuse leidmine, mida on vaja mudelpildi arvutamiseks=========
		!
		do i=1,all_comp%N_comp
			all_comp%comp(i)%mass_abs_tol = leia_massi_abs_tol(all_comp%comp(i), images)
		end do
		!
		! ========= komponentide mudelpiltide arvutamised ==============
		!
		if(kas_koik_pildid_samast_vaatlusest) then
			do i=1,size(mudelid, 1)
				if(mudelid(i)%recalc_image) then
		    		!reaalselt vaja yhe korra ainult teha (koord arvutused sisuslielt)...seega mitteoptimaalsus siin  
					call fill_comp_image(all_comp, i, mudelid(i), via_adaptive_im, kas_los, mida_arvutatakse)
					!kui massid fitib teistest eraldi, siis salvestab massi pildid eraldi
					if(kas_fitib_massid_eraldi) then
						do j=1,size(images); 
							if(kas_rakendab_psf) then
								call rakenda_psf(mudelid(i)%mx, images(j)%psf, pilt_psf)
								to_massfit(j)%M(i,:,:) = pilt_psf
							else
								to_massfit(j)%M(i,:,:) = mudelid(i)%mx; 
							end if
							
						end do
					end if
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
		allocate(weights(1:size(images,1),1:all_comp%N_comp)); weights = -1.234 !ehk algselt koigile mingi tobe v22rtus, et hiljem saaks vigasust kontrollida
		allocate(weights_for_single_im(1:all_comp%N_comp))			
		do i=1,size(images)			
			!
			! ================== kaalude arvutamised ========= 
			!
			do j=1,all_comp%N_comp
					weights(i,j) = images(i)%filter%mass_to_obs_unit(all_comp%comp(j)%dist, all_comp%comp(j)%population_name)
			end do
			if(kas_fitib_massid_eraldi) to_massfit(i)%w(:) = weights(i,:) !hilisemaks kasutamiseks kui otsustab eraldi fittida
		end do
		if(any(weights == -1.234)) then
			print*, "Vale M/L sisend"; stop
		end if

		!t2psustus vastavalt sellele, kas fitib massid eraldi	
		if(kas_fitib_massid_eraldi) then
			!lisab kaaludele veel massi kordaja sisemisest fittimisest
			lisakaalud_massile = fiti_massi_kordajad(to_massfit)
			do j=1,all_comp%N_comp
				weights(:,j) = weights(:,j) * lisakaalud_massile(j) !lisab kaaludesse, et peaks v2hem arvutama
			end do
		else
			if(allocated(lisakaalud_massile)) deallocate(lisakaalud_massile)
			allocate(lisakaalud_massile(1:all_comp%N_comp))
			lisakaalud_massile = 1.0 !ehk v2ljund identsusteisendus
		end if

		!
		! ========= mudelpiltide kokkupanek ===========
		!	
		do i=1,size(images)	
			if(allocated(pilt)) deallocate(pilt)
			if(kas_koik_pildid_samast_vaatlusest) then
				weights_for_single_im = weights(i,:)
				pilt = combine_comp_images_to_make_image(mudelid, weights_for_single_im)
			else
				print*, "Not yet implemented in likelihood"
				stop
			end if


			
			if(.false.) then
				if(allocated(pilt_psf)) deallocate(pilt_psf)
				call convolve(pilt, images(i)%psf, pilt_psf)
				print*, "konvoleeritud"
			else
				if(allocated(pilt_psf)) deallocate(pilt_psf)
				allocate(pilt_psf(1:size(pilt,1), 1:size(pilt,2)))
				pilt_psf = pilt
			end if
			
			!output horedamalt
			if(mod(LL_counter,10)==0)then
				call write_matrix_to_fits(pilt_psf, images(i)%output_mdl_file)
				pilt = (images(i)%obs-pilt_psf)
				where (.not.images(i)%mask) pilt = 0.0
				call write_matrix_to_fits(pilt, images(i)%output_diff_file)
				pilt = abs((images(i)%obs-pilt_psf)/images(i)%sigma) !muutuja yle kasutamine
				where (.not.images(i)%mask) pilt = 0.0
				call write_matrix_to_fits(pilt, images(i)%output_rel_diff_file)
			end if

			
			!
			! ======== loglike ise ========
			!
			
			res = res + sum(-1.0*( (pilt_psf-images(i)%obs)**2*0.5/((images(i)%sigma)**2  + (images(i)%sky_noise**2 + abs(images(i)%obs)))), images(i)%mask) 
		end do
		print*, LL_counter, "LL = ", res
! 		print*, ""
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
	function fiti_massi_kordajad(to_massfit) result(res)
		implicit none
		type(masside_arvutamise_tyyp), dimension(:), allocatable :: to_massfit
		integer :: i,j,k,m, iter
		integer :: N_k, N_i, max_iter
		real(rk), dimension(:), allocatable :: res
		real(rk), dimension(:,:), allocatable :: L_km, L0_km, B_km, inv_L_km !kogu chisq, piltidest chisq, barrier osa chisq-st, poordmaatriks
		real(rk), dimension(:), allocatable :: L_k, L0_k, B_k, nihe
		real(rk), dimension(:,:), allocatable :: tmp_pilt
		real(rk) :: lambda, gamma
! 		real(rk) :: test1
! 		integer :: testi

		lambda = 50.0 !m22rab kui t2pselt ei tohi massid nulli minna... voib olla problemaatiline kui on suured hypped iteratsioonide vahel.. 
		gamma = 0.7 !ehk kui kiiresti liigub iteratsioonide vahel
		max_iter = 50
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
		allocate(res(1:N_k)); res = 1.0; 
		if(.not.allocated(massi_kordajad_eelmine)) then
			allocate(massi_kordajad_eelmine(1:N_k)); !algne masside kordajad on 1.0... seal hakkab edasi roomama
		else
			res = massi_kordajad_eelmine !ehk votab algl2hendi eelimise fittimise tulemusest
		end if


! do testi = 1,10
! do iter=1,5; call random_number(res(iter)); res = res; end do
	
		do iter = 1, max_iter
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
			end do

			!
			! ================= edasi liikumine miinimumi poole =================
			!
			inv_L_km = fun_inv_mx(L_km) !poordmaatriksi leidmine, et edasi liikuda
			nihe = 0.0
			massi_kordajad_eelmine = res !salvestab viimase, et testida koonduvust
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
			if(maxval(abs(res-massi_kordajad_eelmine)/res)<massif_fiti_rel_t2psus) then
! 				print*, "t2psus, iter:", minval(abs(res-massi_kordajad_eelmine)/res), iter
				exit
			end if
		end do

! print"(5F,I)", res, iter
! end do
! stop
! 		!test
! 		test1 = 0.0
! 		do i=1,size(to_massfit, 1)
! 			tmp_pilt = 0.0
! 			do j=1,size(to_massfit(i)%w,1)
! 				tmp_pilt = tmp_pilt + to_massfit(i)%w(j)*res(j)*to_massfit(i)%M(j,:,:)
! 			end do
! 			test1 = test1 - sum( (to_massfit(i)%I-tmp_pilt)**2*0.5*to_massfit(i)%inv_sigma2, to_massfit(i)%mask)
! 		end do
! ! 		print*, "test", test1
	end function fiti_massi_kordajad
	
end module likelihood_module