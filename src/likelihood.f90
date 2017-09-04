module likelihood_module
	use images_module
	use fill_model_image_module
	use file_operations_module !vaja ainult ajutiselt testimisel
	use konvolutsioon_module
	
contains
	function calc_log_likelihood(all_comp, images) result(res)
		implicit none
		type(all_comp_type), intent(inout) :: all_comp
		type(image_type), dimension(:), allocatable, intent(in) :: images
		type(model_image_real_type), dimension(:), allocatable :: mudelid
		real(rk) :: res
		real(rk), dimension(:,:), allocatable :: pilt, pilt_psf
		integer :: i, j, k
		logical :: kas_koik_pildid_samast_vaatlusest !ehk kui koik sama cutout (ruumiliselt ja lahutuselt), siis mudelpilti peab v2he ymber arvutama
		logical :: via_comp_im, kas_los !pildi t2itmise erip2rad
		character(len=default_character_length) :: mida_arvutatakse
		real(rk), dimension(:), allocatable :: weights !ehk M/L suhted fotomeetria korral
		real(rk) :: yhikute_kordaja !10e10Lsun to counts/s
		
		!need peaks tulema mujalt seadetest, mitte k2sitsi
		kas_koik_pildid_samast_vaatlusest = .true. 
		via_comp_im = .true.
		kas_los = .true.
		mida_arvutatakse = "Not in use"
		
! 		print*, "calc log likelihood", all_comp%comp(:)%comp_image_number
		!
		! ========= komponentide mudelpiltide arvutamised ==============
		!
		if(kas_koik_pildid_samast_vaatlusest) then
			allocate(mudelid(1:all_comp%N_comp))
			do i=1,size(mudelid, 1)
				call create_model_image_from_obs(mudelid(i), images(1)) !reaalselt vaja yhe korra ainult teha (koord arvutused sisuslielt)...seega mitteoptimaalsus siin
				call fill_model_image(all_comp, i, mudelid(i), via_comp_im, kas_los, mida_arvutatakse)
! 				print*, all_comp%comp(1)%incl ,sum(mudelid(i)%mx)
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
				do k=1,size(images(i)%filter%population_names, 1)
					if(trim(images(i)%filter%population_names(k)) == trim(all_comp%comp(j)%population_name)) then
! 						weights(j) = images(i)%filter%population_mass_to_light_ratios(k)
						!
						!TODO yhikute kordaja peaks solutma komponendi kaugusets, mis igal comp eraldi
						!
						yhikute_kordaja = 10**(0.4*(images(i)%filter%ZP-images(i)%filter%Mag_sun) + 6.0 - 2.0*log10( all_comp%comp(1)%dist ) )
						weights(j) = yhikute_kordaja / images(i)%filter%population_mass_to_light_ratios(k) !eesm2rk saada pildi yhikutesse korrutades
						exit
					end if
				end do
			end do
			if(any(weights == -1.234)) then
				print*, "Vale M/L sisend"; stop
			end if
			!
			! ========= mudelpiltide kokkupanek ===========
			!
			if(allocated(pilt)) deallocate(pilt)
			if(kas_koik_pildid_samast_vaatlusest) then
				pilt = combine_model_images_to_make_image(mudelid, weights)
			else
				print*, "Not yet implemented"
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
			call write_matrix_to_fits(pilt_psf, images(i)%output_mdl_file)
			!
			! ======== loglike ise ========
			!
			
			res = res + sum(-1.0*( (pilt_psf-images(i)%obs)**2*0.5/((images(i)%sigma)**2  + 0.001 + abs(images(i)%obs))), images(i)%mask) 
			
		end do
		
		

		
		
! 		print*, "LL = ", res
	end function calc_log_likelihood
end module likelihood_module