module yldine_matemaatika_module
	use constants_module
	use distributions_module
	public
	interface coordinate_rotation
		module procedure coordinate_rotation_alpha
		module procedure coordinate_rotation_sincos_alpha
	end interface ! coordinate_rotation
	interface interpolate_1D
		module procedure interpolate_1D_arr
		module procedure interpolate_1D_nr
	end interface
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
	elemental function interpolate_1D_arr(distr, x) result(res)
		implicit none
		type(distribution_1D_type), intent(in) :: distr
		real(rk), intent(in) :: x
		integer :: i
		real(rk) :: res
		i = leia_l2him( x)
		!praegu ykskoik kas interpoleerib voi ekstrapoleerib... v2hem t2pne, aga aeglaselt muutuvatele funktsioonidele sobib... kiirem implemnteerida
		if(i>1) then
			res = interpolate_1D(distr%x(i-1), distr%x(i), distr%f(i-1), distr%f(i), x)
		else
			res = interpolate_1D(distr%x(i), distr%x(i+1), distr%f(i), distr%f(i+1), x)
		end if
	contains
		elemental function leia_l2him(x) result(i)
		!annab l2hima indeksi x-le
			implicit none
			real(rk), intent(in) :: x
			integer :: i
			! ======= meetod 1 ===== brute force - aeglane, aga kindel
			i = minloc(abs(distr%x	-x), 1)
			! ======= meetod 2 ====== not implemented
		end function leia_l2him
	end function interpolate_1D_arr
	elemental function interpolate_1D_nr(x0, x1, val0, val1, x) result(res)
		implicit none
		real(rk), intent(in) :: x0, x1, val0, val1, x
		real(rk) :: res
		real(rk) :: tuletis
		tuletis = (val1 - val0)/(x1-x0)
		res = val0 + tuletis*(x-x0)
	end function interpolate_1D_nr

	elemental function Gauss(x, mean, sigma2) result(res)
		implicit none
		real(rk), intent(in) :: x, mean, sigma2
		real(rk) :: res
		res = exp(-0.5*(x-mean)**2/sigma2)/sqrt(2*pi*sigma2)
	end function Gauss
	function fun_inv_mx(A) result(Ainv)
		!poordmaatriksi leidmine voetud http://fortranwiki.org/fortran/show/Matrix+inversion
		! Returns the inverse of a matrix calculated by finding the LU
		! decomposition.  Depends on LAPACK.
		integer, parameter :: dp=kind(1.0d0)
	  real(rk), dimension(:,:), intent(in) :: A
	  real(rk), dimension(size(A,1),size(A,2)) :: Ainv

	  real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
	  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
	  integer :: n, info

	  ! External procedures defined in LAPACK
	  external DGETRF
	  external DGETRI

	  ! Store A in Ainv to prevent it from being overwritten by LAPACK
	  Ainv = A
	  n = size(A,1)

	  ! DGETRF computes an LU factorization of a general M-by-N matrix A
	  ! using partial pivoting with row interchanges.
	  call DGETRF(n, n, Ainv, n, ipiv, info)

	  if (info /= 0) then
	     stop 'Matrix is numerically singular!'
	  end if

	  ! DGETRI computes the inverse of a matrix using the LU factorization
	  ! computed by DGETRF.
	  call DGETRI(n, Ainv, n, ipiv, work, n, info)

	  if (info /= 0) then
	     stop 'Matrix inversion failed!'
	  end if
	end function fun_inv_mx
end module yldine_matemaatika_module