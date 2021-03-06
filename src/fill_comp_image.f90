module fill_comp_image_module
	!use images
	use comp_image_module
	use all_comp_module
	use adaptive_image_real_module
	use los_real_integration_module
	use approx_Einasto_module
	integer, save :: countersees=0, counterv2ljas=0

contains

	subroutine fill_comp_image(all_comp, comp_nr, mdl)
		implicit none
		logical :: via_ci !ehk via_adaptive_im
		type(all_comp_type), intent(inout) :: all_comp !out osa ainult adaptive_im numbri jaoks
		integer :: comp_nr !millist komponenti arvutatakse
		type(comp_image_real_type), intent(inout) :: mdl
		logical :: kas_Einasto_approx
		real(rk) :: Einasto_q_par
	! 	real(rk) :: test1, test2
		interface
			function pind(Xc, Yc) result(res)
				import rk
				real(rk), intent(in) :: Xc, Yc
				real(rk) :: res
			end function pind
			function ruum(R,z,theta) result(res)
				import rk
				real(rk), intent(in) :: R,z
				real(rk), intent(in), optional :: theta
				real(rk) :: res
			end function ruum
		end interface
		procedure(pind), pointer :: f_ptr, f_ptr_otse
		procedure(ruum), pointer :: ruum_ptr
		!comp im asjad
		type(adaptive_image_type), pointer :: adaptive_im	
	! 	integer :: image_number !comp image number
		!muud vidinad
		integer :: i, j

		!
		! ================== kontrollib ja arvutab komponendi koordinaadid uuesti kui vaja =============
		!
		do j=1,size(mdl%mx, 2)
		do i=1,size(mdl%mx, 1)
			!component coodinates
			call all_comp%comp(comp_nr)%XpYp_to_XcYc(mdl%pix(i,j)%Xp_nurgad,mdl%pix(i,j)%Yp_nurgad, mdl%pix(i,j)%Xc_nurgad, mdl%pix(i,j)%Yc_nurgad) !koordinaatide teisendus
			mdl%pix(i,j)%dXc2 =  sqrt((mdl%pix(i,j)%Xc_nurgad(2) - mdl%pix(i,j)%Xc_nurgad(1))**2 +  (mdl%pix(i,j)%Yc_nurgad(2) - mdl%pix(i,j)%Yc_nurgad(1))**2 ) &
							   * sqrt((mdl%pix(i,j)%Xc_nurgad(4) - mdl%pix(i,j)%Xc_nurgad(1))**2 +  (mdl%pix(i,j)%Yc_nurgad(4) - mdl%pix(i,j)%Yc_nurgad(1))**2 )
		end do
		end do

		!
		! ============= valib, kas kasutab kiirendamiseks adaptive_im arvutust ====================
		!

		via_ci = via_adaptive_im
		via_ci = via_ci .and. all_comp%comp(comp_nr)%prof_den%kas_3D !ehk kui on 2D prof, siis ei kasuta kunagi adaptiivset pilti.
		kas_Einasto_approx = via_ci .and. trim(all_comp%comp(comp_nr)%comp_prof_name)=="Einasto" !ehk tolmu ei kasuta
	
		!valitakse, milline funktsioon arvutamiseks
		ruum_ptr => tihedus 
	
		!
		! =========== konstrueeritakse piksli t2itmiseks sobilik funktsioon vastavalt sellele
		! kas tahab voi mitte yle vaatejoone integreerida  =====================
		!
		if(all_comp%comp(comp_nr)%prof_den%kas_3D) then
			f_ptr => los
			f_ptr_otse => los !valib ikkai teise otse, et siseosad oleks t2psemad
		else
			f_ptr => surface
			f_ptr_otse => surface !valib ikkai teise otse, et siseosad oleks t2psemad
		end if
	
		!  lisamehanismid kui on soov arvutada adaptive_im kaudu
		if(via_ci) then
			if(kas_Einasto_approx) then
! 				print*, "kas Einasto approx", kas_Einasto_approx
				call all_comp%comp(comp_nr)%prof_den%get_val("q", Einasto_q_par)
				call fill_approx_Einasto(f_ptr, all_comp%comp(comp_nr)%incl, Einasto_q_par)
				f_ptr => approx_1D_Einasto
			else
				call fill_adaptive_image_real(f_ptr, all_comp%comp(comp_nr)%adaptive_image_number, all_comp%comp(comp_nr)%mass_abs_tol)
				call get_pointer_to_adaptive_image_number_X(adaptive_im, all_comp%comp(comp_nr)%adaptive_image_number)
				f_ptr => vota_adaptive_im_pilt
			end if

		end if

	
		call fill_corners(mdl, f_ptr_otse, f_ptr); ! print*, "Exact"
! 		call fill_corners(mdl, f_ptr_otse, f_ptr); !print*, "optimaalne"!
	
	
		!
		! ========= pikslite t2itmine =========
		!
		do j= 1,size(mdl%mx, 2)
		do i= 1,size(mdl%mx, 1)
			if(kas_kasutab_l2hendit_piksli_v22rtuse_jaoks(mdl%pix(i,j))) then
				mdl%mx(i,j) = mdl%pix(i,j)%get_val(f_ptr, all_comp%comp(comp_nr)%mass_abs_tol)
			else
				mdl%mx(i,j) = mdl%pix(i,j)%get_val(f_ptr_otse, all_comp%comp(comp_nr)%mass_abs_tol)
			end if
		end do
		end do

		if(associated(f_ptr)) nullify(f_ptr)
		if(associated(ruum_ptr)) nullify(ruum_ptr)
		if(kas_Einasto_approx) call delete_approx_Einasto() !m2lu nullimiseks
	contains
		function approx_1D_Einasto(Xc, Yc) result(res)
			implicit none
			real(rk), intent(in) :: Xc, Yc
			real(rk) :: res
			real(rk) :: A
			res = otsi_approx_Einasto(Xc, Yc)
		end function approx_1D_Einasto
		
		function vota_adaptive_im_pilt(Xc, Yc) result(res)
			implicit none
			real(rk), intent(in) :: Xc, Yc
			real(rk) :: res
			res = adaptive_im%get_val(abs(Xc), abs(Yc))
		end function vota_adaptive_im_pilt
	
		function surface(Xc, Yc) result(res)
			implicit none
			real(rk), intent(in) :: Xc, Yc
			real(rk) :: res
			real(rk) :: R, theta
			R = sqrt((Yc*all_comp%comp(comp_nr)%sec_incl)**2 + Xc * Xc)
			theta = atan2(Yc,(Xc*all_comp%comp(comp_nr)%cos_incl)) + all_comp%comp(comp_nr)%theta0
			res = ruum_ptr(R, 0.0_rk, theta)
		end function surface
	
		function los(Xc, Yc) result(res)
			implicit none
			real(rk), intent(in) :: Xc, Yc
			real(rk) :: res
			res = integreerimine_los(ruum_ptr, all_comp%comp(comp_nr), Xc, Yc, all_comp%comp(comp_nr)%prof_den%fun_los_lopmatus(all_comp%comp(comp_nr)%incl)) !comp on vaja sisse anda koordinaatide teisendamiseks
		end function los
	
		function tihedus(R,z,theta) result(res) !koige lihtsam funktsioon esialgu
			implicit none
			real(rk), intent(in) :: R,z
			real(rk), intent(in), optional :: theta
			real(rk) :: res
			res = all_comp%comp(comp_nr)%prof_den%fun_den(R,z,theta)
		end function tihedus

	
	end subroutine fill_comp_image
	
	
	subroutine fill_comp_image_dustplane(all_comp, comp_nr, mdl)
		implicit none
		logical :: via_ci !ehk via_adaptive_im
		type(all_comp_type), intent(inout) :: all_comp !out osa ainult adaptive_im numbri jaoks
		integer :: comp_nr !millist komponenti arvutatakse
		type(comp_image_real_type), intent(inout) :: mdl
		logical :: kas_enne_tasandit
		logical :: kas_Einasto_approx
		real(rk) :: Einasto_q_par

	! 	real(rk) :: test1, test2
		interface
			function pind(Xc, Yc) result(res)
				import rk
				real(rk), intent(in) :: Xc, Yc
				real(rk) :: res
			end function pind
			function ruum(R,z,theta) result(res)
				import rk
				real(rk), intent(in) :: R,z
				real(rk), intent(in), optional :: theta
				real(rk) :: res
			end function ruum
		end interface
		procedure(pind), pointer :: f_ptr, f_ptr_enne, f_ptr_p2rast, f_ptr_otse
		procedure(ruum), pointer :: ruum_ptr
		!comp im asjad
		type(adaptive_image_type), pointer :: adaptive_im_enne, adaptive_im_p2rast	
	! 	integer :: image_number !comp image number
		!muud vidinad
		integer :: i, j

		!
		! ================== kontrollib ja arvutab komponendi koordinaadid uuesti kui vaja =============
		!
		do j=1,size(mdl%mx, 2)
		do i=1,size(mdl%mx, 1)
			!component coodinates
			call all_comp%comp(comp_nr)%XpYp_to_XcYc(mdl%pix(i,j)%Xp_nurgad,mdl%pix(i,j)%Yp_nurgad, mdl%pix(i,j)%Xc_nurgad, mdl%pix(i,j)%Yc_nurgad) !koordinaatide teisendus
! 			mdl%pix(i,j)%dXc2 =  sqrt((mdl%pix(i,j)%Xc_nurgad(2) - mdl%pix(i,j)%Xc_nurgad(1))**2 +  (mdl%pix(i,j)%Yc_nurgad(2) - mdl%pix(i,j)%Yc_nurgad(1))**2 ) &
! 							   sqrt((mdl%pix(i,j)%Xc_nurgad(4) - mdl%pix(i,j)%Xc_nurgad(1))**2 +  (mdl%pix(i,j)%Yc_nurgad(4) - mdl%pix(i,j)%Yc_nurgad(1))**2 )
		end do
		end do
		mdl%pix(:,:)%dXc2 =  sqrt((mdl%pix(1,1)%Xc_nurgad(2) - mdl%pix(1,1)%Xc_nurgad(1))**2 +  (mdl%pix(1,1)%Yc_nurgad(2) - mdl%pix(1,1)%Yc_nurgad(1))**2 ) &
						   * sqrt((mdl%pix(1,1)%Xc_nurgad(4) - mdl%pix(1,1)%Xc_nurgad(1))**2 +  (mdl%pix(1,1)%Yc_nurgad(4) - mdl%pix(1,1)%Yc_nurgad(1))**2 )
		!
		! ============= valib, kas kasutab kiirendamiseks adaptive_im arvutust ====================
		!
		via_ci = via_adaptive_im
		via_ci = via_ci .and. all_comp%comp(comp_nr)%prof_den%kas_3D !ehk kui on 2D prof, siis ei kasuta kunagi adaptiivset pilti.
		kas_Einasto_approx = via_ci .and. trim(all_comp%comp(comp_nr)%comp_prof_name)=="Einasto"
		if(kas_Einasto_approx) then
			call fill_comp_image(all_comp, comp_nr, mdl) !ehk sisuliselt t2idetakse kogu vaatejoon 2ra ning lahutatakse p2rast pool vaatejoont
			mdl%M_p2rast_tasandit = mdl%mx !kohat2itja, et hiljem saaks ymber kasutada muid asju... 
		end if
		
		!valitakse, milline funktsioon arvutamiseks
		ruum_ptr => tihedus 
	
		!
		! =========== konstrueeritakse piksli t2itmiseks sobilik funktsioon vastavalt sellele
		! kas tahab voi mitte yle vaatejoone integreerida  =====================
		!
		if(all_comp%comp(comp_nr)%prof_den%kas_3D) then
			f_ptr => los_dustplane
			f_ptr_otse => los_dustplane !valib ikkagi teise otse, et siseosad oleks t2psemad
		else
			f_ptr => surface
			f_ptr_otse => surface !valib ikkai teise otse, et siseosad oleks t2psemad
		end if
		!  lisamehanismid kui on soov arvutada adaptive_im kaudu
		if(via_ci) then
			if(.not.kas_Einasto_approx) then
				kas_enne_tasandit = .false.
				call fill_adaptive_image_real(f_ptr, all_comp%comp(comp_nr)%adaptive_im_p2rast_tasandit, all_comp%comp(comp_nr)%mass_abs_tol)
				call get_pointer_to_adaptive_image_number_X(adaptive_im_p2rast, all_comp%comp(comp_nr)%adaptive_im_p2rast_tasandit)
			end if
			kas_enne_tasandit = .true.
			call fill_adaptive_image_real(f_ptr, all_comp%comp(comp_nr)%adaptive_im_enne_tasandit, all_comp%comp(comp_nr)%mass_abs_tol)
			call get_pointer_to_adaptive_image_number_X(adaptive_im_enne, all_comp%comp(comp_nr)%adaptive_im_enne_tasandit)
			f_ptr => vota_adaptive_im_pilt
		end if
	
		!
		! ======== siin t2idetakse mudelpildi enne ja p2rast tasandit eraldi... tauga korrutamine ning kokku panemine on hiljem.
		!
		
		kas_enne_tasandit = .true.;
		call fill_corners(mdl, f_ptr_otse, f_ptr);  !ehk kui Einasto l2hend, siis juba nurgad t2idetud 
		do j= 1,size(mdl%mx, 2)
		do i= 1,size(mdl%mx, 1)
			if(kas_kasutab_l2hendit_piksli_v22rtuse_jaoks(mdl%pix(i,j))) then
				mdl%mx(i,j) = mdl%pix(i,j)%get_val(f_ptr, all_comp%comp(comp_nr)%mass_abs_tol) !otse enne massiivi ei saa t2ita, sest nurgad t2idetud adaptive image vale asjaga
			else
				mdl%mx(i,j) = mdl%pix(i,j)%get_val(f_ptr_otse, all_comp%comp(comp_nr)%mass_abs_tol)
			end if
		end do
		end do
		mdl%M_enne_tasandit = mdl%mx 
		mdl%mx = mdl%M_p2rast_tasandit
		!vastavalt sellele, kas on Einasto l2hend saab arvutada p2rast tasandit v22rtused lihtsalt voi keeruliselt
		if(.not.kas_Einasto_approx) then
			kas_enne_tasandit = .false.;
			call fill_corners(mdl, f_ptr_otse, f_ptr);
			do j= 1,size(mdl%mx, 2)
			do i= 1,size(mdl%mx, 1)
				if(kas_kasutab_l2hendit_piksli_v22rtuse_jaoks(mdl%pix(i,j))) then
					mdl%mx(i,j) = mdl%pix(i,j)%get_val(f_ptr, all_comp%comp(comp_nr)%mass_abs_tol)
				else
					mdl%mx(i,j) = mdl%pix(i,j)%get_val(f_ptr_otse, all_comp%comp(comp_nr)%mass_abs_tol)
				end if
			end do
			end do
			mdl%M_p2rast_tasandit = mdl%mx
		else
			mdl%M_p2rast_tasandit = mdl%mx - mdl%M_enne_tasandit !ehk integreerimise rajadega m2ngides saab kiiremini
! 			print*, sum(mdl%mx), sum(mdl%M_p2rast_tasandit), sum(mdl%M_enne_tasandit)
		end if
		
		!m2lu korrastamine igaks juhuks
		if(associated(f_ptr)) nullify(f_ptr)
		if(associated(f_ptr_enne)) nullify(f_ptr_enne)
		if(associated(f_ptr_p2rast)) nullify(f_ptr_p2rast)
		if(associated(ruum_ptr)) nullify(ruum_ptr)
	contains
		
		function vota_adaptive_im_pilt(Xc, Yc) result(res)
			implicit none
			real(rk), intent(in) :: Xc, Yc
			real(rk) :: res
			if(kas_enne_tasandit) then
				res = adaptive_im_enne%get_val(abs(Xc), Yc)
			else
				res = adaptive_im_p2rast%get_val(abs(Xc), Yc)
			end if
		end function vota_adaptive_im_pilt
	
		function surface(Xc, Yc) result(res)
			implicit none
			real(rk), intent(in) :: Xc, Yc
			real(rk) :: res
			real(rk) :: R, theta
			R = sqrt((Yc*all_comp%comp(comp_nr)%sec_incl)**2 + Xc * Xc)
			theta = atan2(Yc,(Xc*all_comp%comp(comp_nr)%cos_incl)) + all_comp%comp(comp_nr)%theta0
			res = ruum_ptr(R, 0.0_rk, theta)
		end function surface
	
		function tihedus(R,z,theta) result(res) !koige lihtsam funktsioon esialgu
			implicit none
			real(rk), intent(in) :: R,z
			real(rk), intent(in), optional :: theta
			real(rk) :: res
			res = all_comp%comp(comp_nr)%prof_den%fun_den(R,z,theta)
		end function tihedus

		function los_dustplane(Xc, Yc) result(res)
			implicit none
			real(rk), intent(in) :: Xc, Yc
			real(rk) :: res
			res = integreerimine_los_dustplane(ruum_ptr, all_comp%comp(comp_nr), Xc, Yc, all_comp%comp(comp_nr)%prof_den%fun_los_lopmatus(all_comp%comp(comp_nr)%incl), kas_enne_tasandit)
		end function los_dustplane
	
	end subroutine fill_comp_image_dustplane
	
	
	function leia_pix_pindala(im, comp) result(res)
		implicit none
		type(image_type), intent(in) :: im
		type(comp_type), intent(in) :: comp
		real(rk) :: xi1, yi1, xi2, yi2, res
		real(rk) :: xp1, yp1, xp2, yp2
		real(rk) :: xc1, yc1, xc2, yc2
		xi1 = 1.0; yi1 = 1.0
		xi2 = 1.0; yi2 = 2.0
		call convert_XiYi_to_XpYp(im, xi1, yi1 , xp1, yp1)
		call convert_XiYi_to_XpYp(im, xi2, yi2 , xp2, yp2)
		call convert_XpYp_to_XcYc(comp, xp1, yp1, xc1, yc1)
		call convert_XpYp_to_XcYc(comp, xp2, yp2, xc2, yc2)
		res = sqrt((xc2-xc1)**2 + (yc2-yc1)**2)
		!kordus vaja seetottu, et pix scale ei pruudi olla vordne eri suundades
		xi1 = 1.0; yi1 = 1.0
		xi2 = 2.0; yi2 = 1.0
		call convert_XiYi_to_XpYp(im, xi1, yi1 , xp1, yp1)
		call convert_XiYi_to_XpYp(im, xi2, yi2 , xp2, yp2)
		call convert_XpYp_to_XcYc(comp, xp1, yp1, xc1, yc1)
		call convert_XpYp_to_XcYc(comp, xp2, yp2, xc2, yc2)
		res = res * sqrt((xc2-xc1)**2 + (yc2-yc1)**2)
	end function leia_pix_pindala

	

end module fill_comp_image_module