module psf_rakendamine_module
	use constants_module
	use file_operations_module
	use konvolutsioon_module
	real(rk), parameter :: psf_sisse_peab_j22ma = 0.99
	procedure(psf_fun), pointer :: rakenda_psf => rakenda_psf_Fourier!rakenda_psf_toores_joud
	interface
		subroutine psf_fun(pilt, psf, res)
			import rk
			implicit none
			real(rk), intent(in), dimension(:,:), allocatable :: pilt
			real(rk), intent(in), dimension(:,:), allocatable :: psf
			real(rk), intent(out), dimension(:,:), allocatable :: res
		end subroutine psf_fun
	end interface
contains
! 	subroutine crop_psf(psf)
! 		implicit none
! 	end subroutine crop_psf
	subroutine rakenda_psf_toores_joud(pilt, psf, res)
		implicit none
		real(rk), intent(in), dimension(:,:), allocatable :: pilt
		real(rk), intent(in), dimension(:,:), allocatable :: psf
		real(rk), intent(out), dimension(:,:), allocatable :: res
		real(rk),  dimension(:,:), allocatable :: tmp
		integer :: i,j
		integer :: N_psf_i, N_psf_j, N_pilt_i, N_pilt_j
		
		N_pilt_i = size(pilt, 1)
		N_pilt_j = size(pilt, 2)
		N_psf_i = size(psf, 1)
		N_psf_j = size(psf, 2)
		if(allocated(res)) deallocate(res)
		allocate(res(1:N_pilt_i, 1:N_pilt_j))
		allocate(tmp(1:(N_pilt_i+N_psf_i), 1:(N_pilt_j+N_psf_j)))						
		res = -1000.0										
		tmp = 0

		!rakendamine, et suurem tulemuspilt
		do i=1,N_psf_i
			do j=1,N_psf_j
				tmp(i:(i+N_pilt_i-1), j:(j+N_pilt_j-1)) = tmp(i:(i+N_pilt_i-1), j:(j+N_pilt_j-1)) + psf(i,j) * pilt
			end do
		end do
		!tagasi oigesse suurusesse
		if(mod(N_psf_i, 2)==1 .and. mod(N_psf_j, 2)==1) then
			i = (N_psf_i-1)/2+1 !keskkoht
			j = (N_psf_j-1)/2+1 !keskkoht
			res = tmp(i:(N_pilt_i+i-1), i:(N_pilt_j+j-1))			
		else
			stop "We assume that psf has odd number of pixels"
		end if		
		
		if(any(isnan(res))) then
			call write_matrix_to_fits(tmp, "na_tmp.fits")
			call write_matrix_to_fits(psf, "na_psf.fits")
			call write_matrix_to_fits(res, "na_res.fits")
			call write_matrix_to_fits(pilt, "na_pilt.fits")
			stop "psf juures tuli NaN"
		end if
		
	end subroutine rakenda_psf_toores_joud
	subroutine rakenda_psf_Fourier(pilt, psf, res)
		implicit none
		real(rk), intent(in), dimension(:,:), allocatable :: pilt
		real(rk), intent(in), dimension(:,:), allocatable :: psf
		real(rk), intent(inout), dimension(:,:), allocatable :: res
		call convolve(pilt, psf, res)
	end subroutine rakenda_psf_Fourier
end module psf_rakendamine_module