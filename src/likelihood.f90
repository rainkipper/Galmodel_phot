module likelihood_module
	use images_module
	use fill_comp_image_module
! 	use file_operations_module !vaja ainult ajutiselt testimisel
! 	use psf_rakendamine_module
	use Newton_Rhapson_module
	use mass_to_lum_image_module
! 	use mynr_amoeba_module
! 	use filters_module
	
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
	real(rk), dimension(:,:), allocatable :: last_ML, ML_vead
	real(rk), dimension(:), allocatable :: last_comp_M
	type(comp_image_real_type), dimension(:), allocatable, private :: mudelid  !siin hoitakse mudeli andmeid ... init_log_likelihood juures pannakse paika
	type(masside_arvutamise_tyyp), dimension(:), allocatable, private :: to_massfit !lihtsustav muutuja
	real(rk), dimension(:), allocatable, private :: massi_kordajad_eelmine !massi kordajad... globaalne muutuja, et j2rgmine loglike arvutamine oleks hea algl2hend votta
	integer, save :: LL_counter = 0 !lihtsalt, mitu LL juba arvutatud
	real(rk) :: alguse_aeg
	interface
		function eri_likelihoodide_kuju(all_comp, images, output_images) result(res)
			import all_comp_type
			import image_type
			import rk
			implicit none
			type(all_comp_type), intent(inout) :: all_comp
			type(image_type), dimension(:), allocatable, intent(in) :: images
			real(rk), dimension(:,:,:), allocatable, intent(out), optional :: output_images !optional output
			real(rk) :: res
		end function eri_likelihoodide_kuju
	end interface
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
				if(.not.allocated(to_massfit(i)%w)) allocate(to_massfit(i)%w(1:all_comp%N_comp))
				if(.not.allocated(to_massfit(i)%I)) allocate(to_massfit(i)%I(1:size(images(i)%obs,1), size(images(i)%obs,2)))
				if(.not.allocated(to_massfit(i)%mask)) allocate(to_massfit(i)%mask(1:size(images(i)%obs,1), size(images(i)%obs,2)))
				if(.not.allocated(to_massfit(i)%M)) allocate(to_massfit(i)%M(1:all_comp%N_comp, 1:size(images(i)%obs,1), size(images(i)%obs,2)))
				if(.not.allocated(to_massfit(i)%inv_sigma2)) allocate(to_massfit(i)%inv_sigma2(1:size(images(i)%obs,1), size(images(i)%obs,2)))
				to_massfit(i)%I = images(i)%obs
				to_massfit(i)%mask = images(i)%mask
				to_massfit(i)%inv_sigma2 = 1.0/images(i)%sigma**2
			end do
		else
			stop "not implemented in init log likelihood"
		end if
		if(mis_fittimise_tyyp == 2) allocate(ML_vead(1:size(images), 1:all_comp%N_comp))
		if(mis_fittimise_tyyp == 1 .or. mis_fittimise_tyyp == 3) allocate(last_comp_M(1:all_comp%N_comp))
		call cpu_time(alguse_aeg) !lihtsalt arvestamiseks kui palju votab keskmiselt aega LL arvutamine
		mudelid(:)%recalc_image = .true.
	end subroutine init_calc_log_likelihood
	function calc_log_likelihood(all_comp, images, output_images) result(res)
		implicit none
		type(all_comp_type), intent(inout) :: all_comp
		type(image_type), dimension(:), allocatable, intent(in) :: images
		real(rk) :: res
		real(rk) :: dt
! 		real(rk), dimension(:), allocatable :: lisakaalud_massile !vaja kui populatsioone fitib
		real(rk), dimension(:,:,:), allocatable, intent(out), optional :: output_images !optional output
		integer, dimension(1:2) :: for_case
		procedure(eri_likelihoodide_kuju), pointer :: f_LL
		LL_counter = LL_counter + 1
		select case(mis_fittimise_tyyp)
			case(1); f_LL => calc_log_likelihood_populations
			case(2); f_LL => calc_log_likelihood_components
			case(3); f_LL => calc_log_likelihood_populations_with_dust
		case default
			stop "Niisugust fittimise tyypi (popul voi ML) pole olemas"
		end select 
		if(present(output_images)) then; res = f_LL(all_comp, images, output_images)
		else; res = f_LL(all_comp, images) ;end if
		call cpu_time(dt)
		if(.not.kas_vaikselt) print*, LL_counter, "LL = ", res, "dt = ",(dt-alguse_aeg)/LL_counter
		if(isnan(res)) stop "err: LL ei tohi olla nan"
	end function calc_log_likelihood
	function calc_log_likelihood_components(all_comp, images, output_images) result(res)
		implicit none
		type(all_comp_type), intent(inout) :: all_comp
		type(image_type), dimension(:), allocatable, intent(in) :: images
		real(rk), dimension(:,:), allocatable :: from_mass_to_lum, ML_kordajad!esimene indeks pilt, teine komponent
		real(rk), dimension(:,:), allocatable :: mudelpilt
		real(rk), dimension(:), allocatable :: algv22rtused, tmp
		real(rk), dimension(:,:), allocatable :: amoeba_algl2hend
		real(rk) :: res
		real(rk), dimension(:,:,:), allocatable, optional  :: output_images
		real(rk), dimension(:,:), allocatable :: tmp_hess
		integer :: i, j,mis_pilt, iter, ii
		!
		! ========== t2psuse leidmine, mida on vaja mudelpildi arvutamiseks=========
		!
		
		do i=1,all_comp%N_comp
			all_comp%comp(i)%mass_abs_tol = leia_massi_abs_tol(all_comp%comp(i), images) * abs_tol_kordaja  !ehk eeldab, et teatud protsent teatud min ML korral
			all_comp%comp(i)%massi_abs_tol_los = all_comp%comp(i)%mass_abs_tol / leia_pix_pindala(images(1), all_comp%comp(1)) 
		end do
		!
		! =========== mudelpiltide arvutamine ==========
		!
! 		if(mudelid(i)%recalc_image) then
		if(kas_koik_pildid_samast_vaatlusest) then
			if(present(output_images)) then !ainult v2ljundi tegemiseks lopus
				allocate(output_images(1:size(images,1), 1:size(images(1)%obs,1), 1:size(images(1)%obs,2)))
				if(.not.allocated(tmp_hess)) allocate(tmp_hess(1:size(images), 1:all_comp%N_comp)) !allokeerimine, et saada pilte paika
				output_images = 0.0 
			end if
			do j=1,size(mudelid, 1)
				call fill_comp_image(all_comp, j, mudelid(j))
				call write_matrix_to_fits(mudelid(j)%mx, trim(all_comp%comp(j)%comp_name)//".fits")
				do i=1,size(images,1)
					if(kas_rakendab_psf) then
						call rakenda_psf(mudelid(j)%mx, images(i)%psf, mudelpilt)
						to_massfit(i)%M(j,:,:) = mudelpilt
					else
						to_massfit(i)%M(j,:,:) = mudelid(j)%mx
					end if
				end do
			end do
		else; stop "not from same observations"
		end if
! 		else; stop "no recalc?"
! 		end if
		
		!
		! ============= komponentide mudelpiltide kombineerimised ================
		!
		allocate(from_mass_to_lum(1:size(images), 1:all_comp%N_comp))
		do i=1,size(images, 1); do j=1,all_comp%N_comp
			from_mass_to_lum(i,j) = images(i)%filter%mass_to_obs_unit(dist = all_comp%comp(j)%dist)
		end do; end do
		!
		! =============== heleduste arvutamisel sisemine fittimine ================
		!
		allocate(ML_kordajad(1:size(images), 1:all_comp%N_comp)); ML_kordajad = 1.0
			allocate(algv22rtused(1:all_comp%N_comp)); 
			call random_number(algv22rtused); algv22rtused = algv22rtused*5.0 !
			do mis_pilt = 1, size(images)
				!siin arvutatakse LM_kordajaid, ML jaoks poordvaartus vaja votta
				ML_kordajad(mis_pilt,:) = optim_NR(algv22rtused, leia_gradient_ML_jaoks, leia_hessian_ML_jaoks) 
				if(present(output_images)) then
					!leiab massile ka vead siia
 					!https://www.physik.hu-berlin.de/de/gk1504/block-courses/autumn-2010/program_and_talks/Verkerke_part3
					tmp_hess = leia_hessian_ML_jaoks(ML_kordajad(mis_pilt,:))
					do j=1,all_comp%N_comp
						ML_vead(mis_pilt, j) = sqrt(1.0/abs(tmp_hess(j,j))) !kui siin negatiivne, siis on LL miinimum hoopis
					end do
				end if
			end do
		if(present(output_images)) ML_vead = 0.5*ML_vead/ML_kordajad**2 !TODO poordv22rtus pole vist p2ris korrektne... 0.5 sellest et LL defineeritakse 2 korda suuremana
		ML_kordajad = 1.0/ML_kordajad !poordv22rtus reaalsete mass-heledus suhete jaoks
		if(.not.allocated(last_ML)) allocate(last_ML(1:size(ML_kordajad, 1), 1:size(ML_kordajad, 2)))
		last_ML = ML_kordajad !salvestab lopptulemuse jaoks
		!
		! ========== p2ris likelihoodi arvutamine
		!
		res = 0.0
		from_mass_to_lum = from_mass_to_lum / ML_kordajad
		do i=1,size(images); 
			if(allocated(mudelpilt)) deallocate(mudelpilt)
			allocate(mudelpilt(1:size(to_massfit(i)%M, 2), 1:size(to_massfit(i)%M, 3))) !saab v2ltida kui eeldada, et koik pildid samast vaatlusest
			mudelpilt = 0.0
			do j=1,all_comp%N_comp
				mudelpilt = mudelpilt + from_mass_to_lum(i,j)*to_massfit(i)%M(j,:,:)
			end do; 
			if(present(output_images) .and. kas_koik_pildid_samast_vaatlusest) output_images(i,:,:) = mudelpilt !ainult viimase v2ljundi jaoks
			res = res - sum( to_massfit(i)%inv_sigma2*(to_massfit(i)%I - mudelpilt)**2 , to_massfit(i)%mask)
		end do
		res = res * 0.5 !et likelihood, mitte chisq
	contains
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
				tmp = 1.0/images(i)%filter%mass_to_obs_unit(comp%dist) !poordvaartus, kuna vaja massi hajumist saada countside hajumisest
				if(tmp<res) then
					res = tmp
				end if
			end do
		
		end function leia_massi_abs_tol
		function leia_LL_ML_jaoks(x) result(res)
			implicit none
			real(rk), dimension(:), intent(in) :: x
			real(rk) :: res
			real(rk), dimension(:,:), allocatable :: mudelpilt
			integer :: j
			real(rk), dimension(:), allocatable :: tmp
			allocate(mudelpilt(1:size(to_massfit(mis_pilt)%M(:,:,:), 2), 1:size(to_massfit(mis_pilt)%M(:,:,:), 3)))
			mudelpilt = 0.0
			allocate(tmp(1:size(x))); tmp = x; do j=1,size(x); if(x(j)<1.0e-7) tmp(j) = 1.0e-7 ;end do
			do j=1,size(images)
				mudelpilt = mudelpilt + abs(tmp(j))*to_massfit(mis_pilt)%M(j,:,:)*from_mass_to_lum(mis_pilt, j) !abs on amoeba lolluste vastu
			end do
			res = -1.0*sum( (to_massfit(mis_pilt)%I - mudelpilt)**2 * to_massfit(mis_pilt)%inv_sigma2, to_massfit(mis_pilt)%mask) + massi_fiti_lambda*sum(log(abs(tmp)))
			!defineeritud vale m2rgiga, et saaks amoeba abil minimeerida
			res = res * -1.0
		end function leia_LL_ML_jaoks
		function leia_gradient_ML_jaoks(x) result(res)
			implicit none
			real(rk), dimension(:), allocatable, intent(in) :: x
			real(rk), dimension(:), allocatable :: res
			real(rk), dimension(:,:), allocatable :: mudelpilt
			integer :: j
			allocate(mudelpilt(1:size(to_massfit(mis_pilt)%M(:,:,:), 2), 1:size(to_massfit(mis_pilt)%M(:,:,:), 3)))
			mudelpilt = 0.0
			do j=1,size(x)
				mudelpilt = mudelpilt + x(j)*to_massfit(mis_pilt)%M(j,:,:)*from_mass_to_lum(mis_pilt, j)
			end do
			if(allocated(res)) deallocate(res); allocate(res(1:size(x)))
			do j=1,size(x)
				res(j) = sum( to_massfit(mis_pilt)%inv_sigma2 * &
							to_massfit(mis_pilt)%M(j,:,:)*from_mass_to_lum(mis_pilt, j) * &
							(to_massfit(mis_pilt)%I - mudelpilt), to_massfit(mis_pilt)%mask) + massi_fiti_lambda/x(j)
			end do
		end function leia_gradient_ML_jaoks
		function leia_hessian_ML_jaoks(x) result(res)
			implicit none
			real(rk), dimension(:), intent(in) :: x
			real(rk), dimension(:,:), allocatable :: res
			integer :: m,n
			
			if(allocated(res)) deallocate(res)
			allocate(res(1:size(x,1), 1:size(x,1))); res = 0.0 
			do m=1,size(x)
				do n=m,size(x)
					res(m,n) = sum( -1.0*to_massfit(mis_pilt)%M(m,:,:)*from_mass_to_lum(mis_pilt, m)*&
										 to_massfit(mis_pilt)%M(n,:,:)*from_mass_to_lum(mis_pilt, n)*&
										 to_massfit(mis_pilt)%inv_sigma2 , to_massfit(mis_pilt)%mask)
					res(n,m) = res(m,n) !symmeetrilise maatriksi lihtsustus
				end do
				res(m,m) = res(m,m) - massi_fiti_lambda/x(m)**2 !barrier term
			end do
			
		end function leia_hessian_ML_jaoks
	end function calc_log_likelihood_components
	
	
	function calc_log_likelihood_populations(all_comp, images, output_images) result(res)
		implicit none
		type(all_comp_type), intent(inout) :: all_comp
		type(image_type), dimension(:), allocatable, intent(in) :: images
		real(rk) :: res
		real(rk), dimension(:,:), allocatable :: pilt, pilt_psf
		integer :: i, j
		character(len=default_character_length) :: mida_arvutatakse
		real(rk), dimension(:,:), allocatable :: weights !ehk M/L suhted fotomeetria korral ... esimene indeks pilt, teine komponent
		real(rk), dimension(:), allocatable :: weights_for_single_im, lisakaalud_massile
		real(rk), dimension(:,:,:), allocatable, intent(out), optional :: output_images !optional output
		
		!need peaks tulema mujalt seadetest, mitte k2sitsi

		mida_arvutatakse = "Not in use"
		
		!
		! ========== t2psuse leidmine, mida on vaja mudelpildi arvutamiseks=========
		!
		do i=1,all_comp%N_comp
			all_comp%comp(i)%mass_abs_tol = leia_massi_abs_tol(all_comp%comp(i), images) * abs_tol_kordaja
			all_comp%comp(i)%massi_abs_tol_los = all_comp%comp(i)%mass_abs_tol / leia_pix_pindala(images(1), all_comp%comp(1)) 
		end do
		!
		! ========= komponentide mudelpiltide arvutamised ==============
		!
		if(kas_koik_pildid_samast_vaatlusest) then
			if(present(output_images)) then !ainult v2ljundi tegemiseks lopus
				allocate(output_images(1:size(images,1), 1:size(images(1)%obs,1), 1:size(images(1)%obs,2)))
				output_images = 0.0 
			end if
			do i=1,size(mudelid, 1)
				if(mudelid(i)%recalc_image) then
		    		!reaalselt vaja yhe korra ainult teha (koord arvutused sisuslielt)...seega mitteoptimaalsus siin  
					call fill_comp_image(all_comp, i, mudelid(i))
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
		if(any(weights == -1.234)) then; print*, "Vale M/L sisend"; stop; end if

		!t2psustus vastavalt sellele, kas fitib massid eraldi	
		if(kas_fitib_massid_eraldi) then
			!lisab kaaludele veel massi kordaja sisemisest fittimisest
			lisakaalud_massile = fiti_massi_kordajad(to_massfit)
			last_comp_M = lisakaalud_massile
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
			
			!output horedamalt
			if(mod(LL_counter,10)==1)then
				if(allocated(pilt_psf)) deallocate(pilt_psf)
				allocate(pilt_psf(1:size(pilt,1), 1:size(pilt,2)))
				call write_matrix_to_fits(pilt, images(i)%output_mdl_file)
				pilt_psf = (images(i)%obs-pilt)
				where (.not.images(i)%mask) pilt_psf = 0.0
				call write_matrix_to_fits(pilt_psf, images(i)%output_diff_file)
				pilt_psf = abs((images(i)%obs-pilt)/images(i)%sigma) !muutuja yle kasutamine
				where (.not.images(i)%mask) pilt_psf = 0.0
				call write_matrix_to_fits(pilt_psf, images(i)%output_rel_diff_file)
			end if

			
			!
			! ======== loglike ise ========
			!
			if(present(output_images)) output_images(i,:,:) = pilt
			res = res + sum((pilt-images(i)%obs)**2/images(i)%sigma**2, images(i)%mask) 
		end do
		res = res * -0.5
	contains
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
	! 		real(rk) :: lambda, gamma
	! 		real(rk) :: test1
	! 		integer :: testi

		
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
					if(kas_barrier) L_k(k) = L_k(k) - massi_fiti_lambda / res(k) !kui piirab masse seestpoolt
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
					if(kas_barrier) L_km(k,k) = L_km(k,k) + massi_fiti_lambda/res(k)**2 !kui piirab masse seestpoolt... sisaldab ainult diagonaalil olevaid elemente
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
					nihe(k) = nihe(k) * massi_fittimise_hyppe_kordaja
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
	end function calc_log_likelihood_populations
	function calc_log_likelihood_populations_with_dust(all_comp, images, output_images) result(res)
		implicit none
		type(all_comp_type), intent(inout) :: all_comp
		type(image_type), dimension(:), allocatable, intent(in) :: images
		real(rk) :: res
		real(rk), dimension(:,:), allocatable :: pilt, pilt_psf
		integer :: i, j
		character(len=default_character_length) :: mida_arvutatakse
		real(rk), dimension(:,:), allocatable :: weights !ehk M/L suhted fotomeetria korral ... esimene indeks pilt, teine komponent
		real(rk), dimension(:), allocatable :: weights_for_single_im, lisakaalud_massile
		real(rk), dimension(:,:,:), allocatable, intent(out), optional :: output_images !optional output
		real(rk) :: tau
		type(population_type), pointer :: pop
		!need peaks tulema mujalt seadetest, mitte k2sitsi

		mida_arvutatakse = "Not in use"
		
		!
		! ========== t2psuse leidmine, mida on vaja mudelpildi arvutamiseks=========
		!
		do i=1,all_comp%N_comp
			all_comp%comp(i)%mass_abs_tol = leia_massi_abs_tol(all_comp%comp(i), images) * abs_tol_kordaja
			all_comp%comp(i)%massi_abs_tol_los = all_comp%comp(i)%mass_abs_tol / leia_pix_pindala(images(1), all_comp%comp(1)) 
		end do
		!
		! ========= komponentide mudelpiltide arvutamised ==============
		!
		if(kas_koik_pildid_samast_vaatlusest) then
			if(present(output_images)) then !ainult v2ljundi tegemiseks lopus
				allocate(output_images(1:size(images,1), 1:size(images(1)%obs,1), 1:size(images(1)%obs,2)))
				output_images = 0.0 
			end if
			do j=1,size(mudelid, 1)
				if(mudelid(j)%recalc_image) then
		    		!reaalselt vaja yhe korra ainult teha (koord arvutused sisuslielt)...seega mitteoptimaalsus siin  
					!call fill_comp_image(all_comp, i, mudelid(i))
					!kontrollib, kas midagi juba sarnast arvutatud.. kui on siis if-i sees juba asendatakse 2ra
					if(.not.kas_sama_juba_arvutatud(j)) call fill_comp_image_dustplane(all_comp, j, mudelid(j)) 
					!kui massid fitib teistest eraldi, siis salvestab massi pildid eraldi
					if(kas_fitib_massid_eraldi) then
						do i=1,size(images); 
							pop => populations(all_comp%comp(j)%population_number)
							mudelid(j)%mx = make_lum_image(mudelid(j), images(i)%filter, pop, all_comp%dust_comp(1), all_comp%comp(j)%dist)
							if(kas_rakendab_psf) then
								call rakenda_psf(mudelid(j)%mx, images(i)%psf, pilt_psf)
								to_massfit(i)%M(j,:,:) = pilt_psf
							else
								to_massfit(i)%M(j,:,:) = mudelid(j)%mx; 
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
		if(any(weights == -1.234)) then; print*, "Vale M/L sisend"; stop; end if

		!t2psustus vastavalt sellele, kas fitib massid eraldi	
		if(kas_fitib_massid_eraldi) then
			!lisab kaaludele veel massi kordaja sisemisest fittimisest
			lisakaalud_massile = fiti_massi_kordajad(to_massfit)
			last_comp_M = lisakaalud_massile !v2ljundi jaoks (lopliku output jaoks) on see koopia
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
			
			!output horedamalt
			if(mod(LL_counter,10)==1)then
				if(allocated(pilt_psf)) deallocate(pilt_psf)
				allocate(pilt_psf(1:size(pilt,1), 1:size(pilt,2)))
				call write_matrix_to_fits(pilt, images(i)%output_mdl_file)
				pilt_psf = (images(i)%obs-pilt)
				where (.not.images(i)%mask) pilt_psf = 0.0
				call write_matrix_to_fits(pilt_psf, images(i)%output_diff_file)
				pilt_psf = abs((images(i)%obs-pilt)/images(i)%sigma) !muutuja yle kasutamine
				where (.not.images(i)%mask) pilt_psf = 0.0
				call write_matrix_to_fits(pilt_psf, images(i)%output_rel_diff_file)
			end if

			
			!
			! ======== loglike ise ========
			!
			if(present(output_images)) output_images(i,:,:) = pilt
			res = res + sum((pilt-images(i)%obs)**2/images(i)%sigma**2, images(i)%mask) 
		end do
		res = res * -0.5
	contains
		function kas_sama_juba_arvutatud(id) result(kas_vaste_olemas)
			!sisuliselt vaatab id j2rgi, kas samade parameetritega komponent on v2iksema id-ga juba olemas.
			!kui on olemas, siis kopeerib pildid jm oige sisse ning annab vasteks TRUE, muidu FALSE
			implicit none
			integer, intent(in) :: id
			logical :: kas_vaste_olemas
			logical :: kas_see_on_sama
			integer :: i, j
			real(rk) :: par1, par2
			character(len=*), dimension(1:3), parameter :: E_par_nimed = ["a0", "q", "N"] !struktuuri parameetrid... mass ei kuulu siia
			
			!
			! ========= otsib vastet =========
			!
			kas_vaste_olemas = .false. !default vastus, et pole ja koik arvutatakse uuesti
			if(id == 1 .or. .not.kas_koik_pildid_samast_vaatlusest) return
			do i=1,id-1
				kas_see_on_sama = .true.
				kas_see_on_sama = kas_see_on_sama .and. (abs(all_comp%comp(i)%dist - all_comp%comp(id)%dist)<epsilon(1.0_rk))
				kas_see_on_sama = kas_see_on_sama .and. (abs(all_comp%comp(i)%pos - all_comp%comp(id)%pos)<epsilon(1.0_rk))
				kas_see_on_sama = kas_see_on_sama .and. (abs(all_comp%comp(i)%theta0 - all_comp%comp(id)%theta0)<epsilon(1.0_rk))
				kas_see_on_sama = kas_see_on_sama .and. (abs(all_comp%comp(i)%incl - all_comp%comp(id)%incl)<epsilon(1.0_rk))
				kas_see_on_sama = kas_see_on_sama .and. (trim(all_comp%comp(i)%comp_prof_name) == trim(all_comp%comp(id)%comp_prof_name))
				if(.not.kas_see_on_sama) cycle !ehk parameetreid ei kontrolli kui midagi juba nihkes
				select case(trim(all_comp%comp(i)%comp_prof_name))
				case("Einasto")
					do j=1,size(E_par_nimed, 1)
						call all_comp%comp(i)%prof_den%get_val(E_par_nimed(j), par1)
						call all_comp%comp(id)%prof_den%get_val(E_par_nimed(j), par2)
						kas_see_on_sama = kas_see_on_sama .and. abs(par1-par2)<epsilon(1.0_rk)
					end do
				case default
					print*, "Kiirendavat vordsust ei ole, arvutab asju palju uuesti (voibolla)"
					kas_see_on_sama = .false.
				end select
				if(kas_see_on_sama) then
					kas_vaste_olemas = .true.
					exit
				else
					return !ehk kui vastet pole, siis v2lja subroutine-ist
				end if
			end do
			!
			! ========== kui vaste olemas, siis kopeerib vajalikud asjad ============
			!
			 mudelid(id)%M_p2rast_tasandit = mudelid(i)%M_p2rast_tasandit 
			 mudelid(id)%M_enne_tasandit = mudelid(i)%M_enne_tasandit 
			 all_comp%comp(id)%adaptive_im_enne_tasandit = all_comp%comp(i)%adaptive_im_enne_tasandit
			 all_comp%comp(id)%adaptive_im_p2rast_tasandit = all_comp%comp(i)%adaptive_im_p2rast_tasandit
			 !nurkade koordinaatide infot pole vaja arvutada kuna pikslid juba arvutatud.
			 if(.false.) then
				 mudelid(id)%pix(:,:)%dXc2 = 0.0 !muud koordinaadid oleka ka muidu vaja
			 end if
		end function kas_sama_juba_arvutatud
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
	! 		real(rk) :: lambda, gamma
	! 		real(rk) :: test1
	! 		integer :: testi

		
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
					if(kas_barrier) L_k(k) = L_k(k) - massi_fiti_lambda / res(k) !kui piirab masse seestpoolt
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
					if(kas_barrier) L_km(k,k) = L_km(k,k) + massi_fiti_lambda/res(k)**2 !kui piirab masse seestpoolt... sisaldab ainult diagonaalil olevaid elemente
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
					nihe(k) = nihe(k) * massi_fittimise_hyppe_kordaja
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
	end function calc_log_likelihood_populations_with_dust

	
	
end module likelihood_module