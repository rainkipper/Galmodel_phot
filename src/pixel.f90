module pixel_module
	use constants_module

	real(rk), parameter :: sqrt2 = sqrt(2.0_rk)
	real(rk), parameter :: edasi_jagamise_rel_t2psus = 0.005
	real(rk), parameter :: edasi_jagamise_abs_t2psus = 0.000
	integer, parameter  :: maxlevel = 5
	type :: square_pixel_type
		real(rk) :: val
		real(rk), dimension(1:4) :: Xi_nurgad
		real(rk), dimension(1:4) :: Yi_nurgad
		real(rk), dimension(1:4) :: Xp_nurgad
		real(rk), dimension(1:4) :: Yp_nurgad
		logical :: recalc_XcYc_coords = .true. !kas vaja koordinaadid yle arvutada
		real(rk), dimension(1:4) :: Xc_nurgad !komponendi koordinaatides, et oleks kiirem arvutada... vb ei ole vaja
		real(rk), dimension(1:4) :: Yc_nurgad
		real(rk), dimension(1:4) :: val_nurgad
		real(rk) :: dXc2 !ehk kylje pikkuse ruut komponendi koordinaatides
	contains
		procedure :: get_val  => get_val_sq
	end type square_pixel_type
	

contains

	
	function get_val_sq(pix, func) result(res)
		class(square_pixel_type), intent(inout) :: pix
! 		type(square_pixel_type) :: pix
		interface 
			function func(Xc, Yc) result(res)
				import rk
				real(rk), intent(in) :: Xc, Yc
				real(rk) :: res
			end function func
		end interface
		real(rk) :: res1, res2, res
		real(rk), dimension(1:3) :: x,y,val
		x(1) = pix%Xc_nurgad(1); x(2) = pix%Xc_nurgad(2); x(3) = pix%Xc_nurgad(4);
		y(1) = pix%Yc_nurgad(1); y(2) = pix%Yc_nurgad(2); y(3) = pix%Yc_nurgad(4);
		val(1) = pix%val_nurgad(1); val(2) = pix%val_nurgad(2); val(3) = pix%val_nurgad(4);
		call t2isnurkse_vordkylgse_kolmnurga_integraal(pix%dXc2, val(1), val(2), val(3), res1)
		call leia_kolmnurga_v22rtus(x,y,val, res1, 1)
		x(1) = pix%Xc_nurgad(3); y(1) = pix%Yc_nurgad(3); val(1) = pix%val_nurgad(3); !2 ylej22nud nurka j22vad samaks
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
			integer :: level !n2itab, mitu iteratsiooni sygaval on... vaja selleks, et lopmatult ei hakkaks jama arvutama
			!
			! ================ alamkolmnurkade integraalide arvutamine =================
			!
			dx2 = (x(2)-x(1))**2 + (y(2)-y(1))**2
			xnext(1) = (x(2)+x(3))*0.5; ynext(1) = (y(2)+y(3))*0.5
			val_t2isnurgas = func(xnext(1), ynext(1)) !tuleb pealmisest subroutine-st
			xnext(3) = x(1); ynext(3) = y(1)
			
			xnext(2) = x(2); ynext(2) = y(2)
			call t2isnurkse_vordkylgse_kolmnurga_integraal(dx2, val_t2isnurgas, val(2), val(1), integraal1)
			xnext(2) = x(3); ynext(2) = y(3)
			call t2isnurkse_vordkylgse_kolmnurga_integraal(dx2, val_t2isnurgas, val(3), val(1), integraal2)
			!
			! =========== vaatab, kas peab edasi arvutama =========== 
			!
			tmp = abs(integraal1 + integraal2 - integraal) !t2psus
			if(level<maxlevel .and. ( (tmp >  edasi_jagamise_abs_t2psus) .or. (abs(tmp/integraal) > edasi_jagamise_rel_t2psus))) then
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
			real(rk), intent(in) :: dx2
			real(rk), intent(out) :: res
			res = (dx2*0.5) / 6.0*(fx+fy-f0)
		end subroutine t2isnurkse_vordkylgse_kolmnurga_integraal
	end function get_val_sq
	
end module pixel_module