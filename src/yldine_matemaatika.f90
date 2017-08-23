module yldine_matemaatika_module
	use constants_module
	public
	interface coordinate_rotation
		module procedure coordinate_rotation_alpha
		module procedure coordinate_rotation_sincos_alpha
	end interface ! coordinate_rotation
contains
	elemental subroutine coordinate_rotation_alpha(xin, yin, alpha, xout, yout)
		implicit none
		real(rk), intent(in) :: xin, yin, alpha
		real(rk), intent(out) :: xout, yout
		real(rk) :: s,c
		s = sin(alpha)
		c = cos(alpha)
		xout = xin*c - yin*s
		yout = xin*s + yin*c
	end subroutine coordinate_rotation_alpha
	elemental subroutine coordinate_rotation_sincos_alpha(xin, yin, sin_alpha, cos_alpha, xout, yout)
		implicit none
		real(rk), intent(in) :: xin, yin, sin_alpha, cos_alpha
		real(rk), intent(out) :: xout, yout
		xout = xin*cos_alpha - yin*sin_alpha
		yout = xin*sin_alpha + yin*cos_alpha
	end subroutine coordinate_rotation_sincos_alpha
	
	elemental function interpolate_1D(x0, x1, val0, val1, x) result(res)
		implicit none
		real(rk), intent(in) :: x0, x1, val0, val1, x
		real(rk) :: res
		real(rk) :: tuletis
		tuletis = (val1 - val0)/(x1-x0)
		res = val0 + tuletis*(x-x0)
	end function interpolate_1D
	elemental function Gauss(x, mean, sigma2) result(res)
		implicit none
		real(rk), intent(in) :: x, mean, sigma2
		real(rk) :: res
		res = exp(-0.5*(x-mean)**2/sigma2)/sqrt(2*pi*sigma2)
	end function Gauss
end module yldine_matemaatika_module