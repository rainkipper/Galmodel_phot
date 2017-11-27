module pixel_module
	use constants_module
	type :: square_pixel_type
		real(rk) :: val
		real(rk), dimension(1:4) :: Xi_nurgad !jrk on 1=all vasak, 2=ylal vasak, 3=ylal_parem, 4=all_parem
		real(rk), dimension(1:4) :: Yi_nurgad
		real(rk), dimension(1:4) :: Xp_nurgad
		real(rk), dimension(1:4) :: Yp_nurgad
		logical :: recalc_XcYc_coords = .true. !kas vaja koordinaadid yle arvutada
		real(rk), dimension(1:4) :: Xc_nurgad !komponendi koordinaatides, et oleks kiirem arvutada... vb ei ole vaja
		real(rk), dimension(1:4) :: Yc_nurgad
		real(rk), dimension(1:4) :: val_nurgad
		real(rk) :: dXc2 !ehk kylje pikkuse ruut komponendi koordinaatides
	contains
! 		procedure :: get_val  => get_val_sq_v2ga_j2me
		procedure :: get_val  => get_val_sq_kolmnurgad
	end type square_pixel_type
	
! 	integer, parameter :: minlevel = 3
contains
	
	elemental function kas_kasutab_l2hendit_piksli_v22rtuse_jaoks(pix) result(res)
		implicit none
		type(square_pixel_type), intent(in) :: pix
		logical :: res
		!vaatab kasti kujuliselt
		res = .not. any( (abs(pix%Xc_nurgad)<adaptive_image_dist_piirang) .and. (abs(pix%Yc_nurgad)<adaptive_image_dist_piirang))
	end function kas_kasutab_l2hendit_piksli_v22rtuse_jaoks

	function get_val_sq_v2ga_j2me(pix, func, abs_tol) result(res)
		implicit none
		class(square_pixel_type), intent(inout) :: pix
		interface 
			function func(Xc, Yc) result(res)
				import rk
				real(rk), intent(in) :: Xc, Yc
				real(rk) :: res
			end function func
		end interface
		real(rk), intent(in) :: abs_tol
		real(rk) :: res
		
		res = sum(pix%val_nurgad)/size(pix%val_nurgad, 1) * pix%dXc2 !ehk keskmine korda pindala
		pix%val = res + 0*abs_tol
	end function get_val_sq_v2ga_j2me
	
	
	
	
	function get_val_sq_kolmnurgad(pix, func, abs_tol) result(res)
		class(square_pixel_type), intent(inout) :: pix
		real(rk), intent(in), optional :: abs_tol
		interface 
			function func(Xc, Yc) result(res)
				import rk
				real(rk), intent(in) :: Xc, Yc
				real(rk) :: res
			end function func
		end interface
		real(rk) :: res1, res2, res
		real(rk) :: edasi_jagamise_abs_t2psus
		real(rk), dimension(1:3) :: x,y,val
		if(present(abs_tol)) then
			edasi_jagamise_abs_t2psus = abs_tol
		else
			edasi_jagamise_abs_t2psus = 1.0e10 !ehk ei arvesta
		end if
		x(1) = pix%Xc_nurgad(1); x(2) = pix%Xc_nurgad(2); x(3) = pix%Xc_nurgad(4);
		y(1) = pix%Yc_nurgad(1); y(2) = pix%Yc_nurgad(2); y(3) = pix%Yc_nurgad(4);
		val(1) = pix%val_nurgad(1); val(2) = pix%val_nurgad(2); val(3) = pix%val_nurgad(4);
		call t2isnurkse_vordkylgse_kolmnurga_integraal(pix%dXc2, val(1), val(2), val(3), res1) !algse l2hendi leidmine
		call leia_kolmnurga_v22rtus(x,y,val, res1, 1) !l2hendi t2psustamine
		x(1) = pix%Xc_nurgad(3); y(1) = pix%Yc_nurgad(3); val(1) = pix%val_nurgad(3); !kaks ylej22nud nurka j22vad samaks
		call t2isnurkse_vordkylgse_kolmnurga_integraal(pix%dXc2, val(1), val(2), val(3), res2)
		call leia_kolmnurga_v22rtus(x,y,val, res2, 1)

		res = res1 + res2
		pix%val = res
		
	contains
		recursive subroutine leia_kolmnurga_v22rtus(x,y,val, integraal, level)
			implicit none
			real(rk), intent(in), dimension(1:3) :: x,y,val !1 on t2isnurga serv, 2 ja 3 ylej22nud
			real(rk), intent(inout) :: integraal
			real(rk), dimension(1:3) :: xnext, ynext, valnext
			real(rk) :: val_t2isnurgas
			real(rk) :: integraal1, integraal2
			real(rk) :: dx2
			real(rk) :: tmp
			logical :: kas_kolm_punkti_sarnased
			logical :: kas_jagab_edasi
			integer :: level !n2itab, mitu iteratsiooni sygaval on... vaja selleks, et lopmatult ei hakkaks jama arvutama
			integer :: loc1 !mediaani leidmiseks
! 			print*, val
! print*, "kolmnurgas"
! 			loc1 = minloc( abs((abs( val - minval(val,1) )/maxval( val - minval(val,1),1 )) - 0.5 ), 1 )
			kas_kolm_punkti_sarnased = sum(val)-maxval(val) > 1.0*maxval(val)
! 			kas_kolm_punkti_sarnased = abs(sum(val)/3.0 - val(loc1)) <  0.2*sqrt( sum((val - sum(val)/3.0)**2)/3.0 )
! 			print*, kas_kolm_punkti_sarnased
! 			if(minval(abs(x),1)<0.2 .and. minval(abs(y),1)<0.2) print*, kas_kolm_punkti_sarnased
! 			print*, loc1
! 			print*, "----"
! 			kas_kolm_punkti_sarnased =
			!
			! ================ alamkolmnurkade integraalide arvutamine =================
			!
			dx2 = ((x(2)-x(1))**2 + (y(2)-y(1))**2)/2.0 !siin jagammine, kuna j2rgmise iteratsiooni pindala vaja
			xnext(1) = (x(2)+x(3))*0.5; ynext(1) = (y(2)+y(3))*0.5
			val_t2isnurgas = func(xnext(1), ynext(1)) !func tuleb pealmisest subroutine-st
			xnext(3) = x(1); ynext(3) = y(1)
			
			xnext(2) = x(2); ynext(2) = y(2)
			call t2isnurkse_vordkylgse_kolmnurga_integraal(dx2, val_t2isnurgas, val(2), val(1), integraal1)
			xnext(2) = x(3); ynext(2) = y(3)
			call t2isnurkse_vordkylgse_kolmnurga_integraal(dx2, val_t2isnurgas, val(3), val(1), integraal2)
			!
			! =========== vaatab, kas peab edasi arvutama =========== 
			!
			tmp = abs(integraal1 + integraal2 - integraal) !t2psus
			
			kas_jagab_edasi = tmp >  edasi_jagamise_abs_t2psus 
			kas_jagab_edasi = kas_jagab_edasi .or. (.not.kas_kolm_punkti_sarnased .and. (tmp/integraal) > pix_edasi_jagamise_rel_t2psus)
			kas_jagab_edasi = kas_jagab_edasi .and. level<pix_iter_maxlevel


			if(  kas_jagab_edasi) then
				valnext(1) = val_t2isnurgas; valnext(3) = val(1);  !molemal juhul samad
				valnext(2) = val(3); xnext(2) = x(3); ynext(2) = y(3)
				call leia_kolmnurga_v22rtus(xnext, ynext, valnext, integraal2, level+1)
				valnext(2) = val(2); xnext(2) = x(2); ynext(2) = y(2)
				call leia_kolmnurga_v22rtus(xnext, ynext, valnext, integraal1, level+1)
			end if
			integraal = integraal1 + integraal2 !summa ikkagi t2psem			
		end subroutine leia_kolmnurga_v22rtus
		subroutine t2isnurkse_vordkylgse_kolmnurga_integraal(dx2, f0, fx, fy, res)
			implicit none
			real(rk), intent(in) :: f0, fx, fy
			real(rk), intent(in) :: dx2!dx2 on siin kylje pikkuse ruut, mitte mingi pindala
			real(rk), intent(out) :: res
			res = (dx2) / 6.0*(fx+fy+f0) 
		end subroutine t2isnurkse_vordkylgse_kolmnurga_integraal
	end function get_val_sq_kolmnurgad
	
end module pixel_module