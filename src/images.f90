module images_module
	use constants_module
	use filters_module
	use yldine_matemaatika_module
	
	type image_type
		!piltide sisend
		character(len=default_character_length)			:: name
		
		type(filter_type), pointer 						:: filter
		character(len=default_character_length)			:: obs_file !sisendfaili nimi
		character(len=default_character_length)			:: mask_file !maskifaili nimi
		character(len=default_character_length)			:: sigma_file !kaalufaili nimi
		character(len=default_character_length)			:: psf_file !psf faili nimi	
		real(rk), dimension(:,:), allocatable 			:: obs
		logical, dimension(:,:), allocatable 			:: mask
		real(rk), dimension(:,:), allocatable 			:: sigma
		real(rk), dimension(:,:), allocatable 			:: psf

		!j2rgnevad peaks olema m22ratud vaatlustest kauguse t2psusega
		real(rk)										:: x !pildi enda nihe mudeli kasutatavate koordinaatide suhtes.... 
		real(rk)										:: y !pildi enda nihe mudeli kasutatavate koordinaatide suhtes. 
		real(rk)										:: pos !pildi enda nihe mudeli kasutatavate koordinaatide suhtes. 
		real(rk) 										:: scale_x !... arcsec = scale * pix
		real(rk) 										:: scale_y !... arcsec = scale * pix
		

		
		
		!v2ljund
		character(len=default_character_length)			:: output_mdl_file
		character(len=default_character_length)			:: output_diff_file
		character(len=default_character_length)			:: output_rel_diff_file
		
		
		!automaatselt m22ratavad muutujad
		real(rk) 										:: cos_pos, sin_pos, tan_pos, sec_pos !vb pole vaja
	contains
		procedure :: XiYi_to_XpYp => convert_XiYi_to_XpYp
	end type image_type
	

contains
	elemental subroutine convert_XiYi_to_XpYp(im, Xi, Yi, Xp, Yp)
		implicit none
		class(image_type), intent(in) :: im
		real(rk), intent(in) :: Xi, Yi
		real(rk), intent(out) :: Xp, Yp
		real(rk) :: x,y
! 		call coordinate_rotation_sincos_alpha(xin, yin, sin_alpha, cos_alpha, xout, yout)
		x = (Xi - im%x) * im%scale_x * arcsec_to_rad
		y = (Yi - im%y) * im%scale_y * arcsec_to_rad
		call coordinate_rotation(x, y, im%sin_pos, im%cos_pos, Xp, Yp)
	end subroutine convert_XiYi_to_XpYp
end module images_module