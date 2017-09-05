module fill_model_image_module
	!use images
	use model_image_module
	use all_comp_module
	use comp_image_real_module
	use los_real_integration_module


contains


subroutine fill_model_image(all_comp, comp_nr, mdl, via_comp_im, kas_los, mida_arvutatakse)
	implicit none
	logical, intent(in) :: kas_los !kas integreerib yle vaatejoone voi votab z=0 v22rtuse
	logical, intent(in), optional :: via_comp_im
	logical :: via_ci !ehk via_comp_im
	type(all_comp_type), intent(inout) :: all_comp !out osa ainult comp_im numbri jaoks
	integer :: comp_nr !millist komponenti arvutatakse
	type(model_image_real_type), intent(inout) :: mdl
	character(len=default_character_length), intent(in), optional :: mida_arvutatakse !ei kasuta... praegu votab alati tiheduses
	real(rk) :: test1, test2
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
	procedure(pind), pointer :: f_ptr
	procedure(ruum), pointer :: ruum_ptr
	!comp im asjad
	type(comp_image_type), pointer :: comp_im	
	integer :: image_number !comp image number
	!muud vidinad
	integer :: i, j
! 	print*, "fill model image", all_comp%comp(:)%comp_image_number
	!
	! ================== kontrollib ja arvutab komponendi koordinaadid uuesti kui vaja =============
	!
	do j=1,size(mdl%mx, 2)
	do i=1,size(mdl%mx, 1)
		!component coodinates
		call all_comp%comp(comp_nr)%XpYp_to_XcYc(mdl%pix(i,j)%Xp_nurgad,mdl%pix(i,j)%Yp_nurgad, mdl%pix(i,j)%Xc_nurgad, mdl%pix(i,j)%Yc_nurgad)
		mdl%pix(i,j)%dXc2 = (mdl%pix(i,j)%Xc_nurgad(3) - mdl%pix(i,j)%Xc_nurgad(1))**2 + (mdl%pix(i,j)%Yc_nurgad(3) - mdl%pix(i,j)%Yc_nurgad(1))**2
	end do
	end do

	!
	! ============= valib, kas kasutab kiirendamiseks comp_im arvutust ====================
	!
	if(.not.present(via_comp_im)) then
		via_ci = .false.
	else
		via_ci = via_comp_im
		!TODO kui comp_im juba olemas, siis siit peaks saama ka millise numbri votab ning kas arvutab ymber
	end if
	
	!valitakse, milline funktsioon arvutamiseks
	ruum_ptr => tihedus 
	
	!
	! =========== konstrueeritakse piksli t2itmiseks sobilik funktsioon vastavalt sellele
	! kas tahab voi mitte yle vaatejoone integreerida  =====================
	!
	if(kas_los) then
		f_ptr => los
	else
		f_ptr => surface
	end if
	
	!  lisamehanismid kui on soov arvutada comp_image kaudu
	if(via_ci) then
		call fill_comp_image_real(f_ptr, all_comp%comp(comp_nr)%comp_image_number, all_comp%comp(comp_nr)%mass_abs_tol)
		call get_pointer_to_comp_im_number_X(comp_im, all_comp%comp(comp_nr)%comp_image_number)
		f_ptr => vota_comp_im_pilt
	end if
	
	call fill_corners_exact(mdl, f_ptr)
! 	call fill_corners(mdl, f_ptr) ... ei toimi praegu
	
	
	!
	! ========= pikslite t2itmine =========
	!
	do j=1,size(mdl%mx, 2)
	do i=1,size(mdl%mx, 1)
		mdl%mx(i,j) = mdl%pix(i,j)%get_val(f_ptr, all_comp%comp(comp_nr)%mass_abs_tol)
	end do
	end do
! 	print*, "fill_model_image", all_comp%comp(1)%sec_incl,all_comp%comp(1)%incl, sum(mdl%mx)
	
contains
	
	function vota_comp_im_pilt(Xc, Yc) result(res)
		implicit none
		real(rk), intent(in) :: Xc, Yc
		real(rk) :: res
		res = comp_im%get_val(Xc, Yc)
	end function vota_comp_im_pilt
	
	function surface(Xc, Yc) result(res)
		implicit none
		real(rk), intent(in) :: Xc, Yc
		real(rk) :: res
		real(rk) :: R
		!
		!TODO kui on mittetelgsymmeetriline, siis theta peaks ka panuse andma... praegu seda ei arvesta
		!
		R = sqrt((Yc*all_comp%comp(comp_nr)%sec_incl)**2 + Xc * Xc)
		res = ruum_ptr(R, 0.0_rk, 0.0_rk)
	end function surface
	
	function los(Xc, Yc) result(res)
		implicit none
		real(rk), intent(in) :: Xc, Yc
		real(rk) :: res
		res = integreerimine_los(ruum_ptr, all_comp%comp(comp_nr), Xc, Yc) !comp on vaja sisse anda koordinaatide teisendamiseks
	end function los
	
	function tihedus(R,z,theta) result(res) !koige lihtsam funktsioon esialgu
		implicit none
		real(rk), intent(in) :: R,z
		real(rk), intent(in), optional :: theta
		real(rk) :: res
		res = all_comp%comp(comp_nr)%prof_den%fun_den(R,z,theta)
	end function tihedus

	
end subroutine fill_model_image


	

end module fill_model_image_module