module psf_rakendamine_module
	use constants_module
	use file_operations_module
	use konvolutsioon_module
	procedure(psf_fun), pointer :: rakenda_psf => rakenda_psf_toores_joud!rakenda_psf_toores_joud
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
	subroutine testi_psf(pilt, psf)
		implicit none
		real(rk), intent(in), dimension(:,:), allocatable :: pilt
		real(rk), intent(in), dimension(:,:), allocatable :: psf
		real(rk), dimension(:,:), allocatable :: tmp
		integer :: i,N
		real(rk) :: t1,t2
! 		call write_matrix_to_fits(pilt, "algne.fits")
! 		call write_matrix_to_fits(psf, "psf.fits")
N = 30
		call cpu_time(t1)
		do i=1,N
			call rakenda_psf_toores_joud(pilt, psf, tmp)
		end do
		call cpu_time(t2)
		print*, "psf suurus", size(psf,1)
		print*, "Toores joud", t2-t1
! 		call write_matrix_to_fits(tmp, "toores.fits")
		do i=1,N
			call rakenda_psf_Fourier(pilt, psf, tmp)
		end do
		call cpu_time(t1)
		print*, "Fourier", t1-t2
! 		call write_matrix_to_fits(tmp, "fourier.fits")
		stop "psf test tehtud"
	end subroutine testi_psf
	subroutine crop_psf(psf, crop, check)
		implicit none
		real(rk), dimension(:,:), allocatable, intent(inout) :: psf
		logical, intent(in), optional :: check, crop
		integer :: cntx, cnty, cnt(1:2)
		integer :: i, N
		real(rk) :: kogulum, siselum
		real(rk), dimension(:,:), allocatable :: pisipsf
		N = size(psf, 1)
		cnt = maxloc(psf); cntx = cnt(1); cnty = cnt(2)
		!check if psf is suitable for modelling assumptions
		if((present(check) .and. check) .or. (.not.present(check))) then
			if(N.ne.size(psf,2) ) stop "psf ei ole ruudu kujuline"
			if(mod(N,2).ne.1) stop "psf ei ole paaritu arv piksleid suur"
			if((cntx .ne. cnty) .or. cntx .ne. ((N-1)/2+1) ) stop "psf heledaim piksel ei ole tsentris"
			!kontoll, kas on heledaim piirkond keskel... ehk keskel peab igas suunas olema langus
			if( 0.5*( psf(cntx-1, cnty) + psf(cntx, cnty) ) < psf(cntx+1, cnty) .or. &
				0.5*( psf(cntx+1, cnty) + psf(cntx, cnty) ) < psf(cntx-1, cnty) .or. &
				0.5*( psf(cntx, cnty-1) + psf(cntx, cnty) ) < psf(cntx, cnty+1) .or. &
				0.5*( psf(cntx, cnty+1) + psf(cntx, cnty) ) < psf(cntx, cnty-1)) then
				stop "psfil ei ole tsentris maksimumi"
			end if
		end if
		!cropping part
		if((present(crop) .and. crop) .or. .not.present(crop)) then
			kogulum = sum(psf)
			do i=1,(N-1)/2+1
				siselum = sum(psf(cntx-i:cntx+i, cnty-i:cnty+i)) !vaatab kui palju j22b teatud raadiuse sisse (ruuduline)
				if(siselum>kogulum*psf_sisse_peab_j22ma) then !psf_sisse_peab_j22ma on defineeritud constants moodulis
					allocate(pisipsf(1:1+2*i, 1:1+2*i))
					pisipsf = psf(cntx-i:cntx+i, cnty-i:cnty+i)
					exit
				end if
			end do
			deallocate(psf)
			allocate(psf(1:size(pisipsf,2), 1:size(pisipsf,2)))
			psf = pisipsf
		end if
		psf = psf / sum(psf) !normeerimine, et massikadu ei oleks... vaja teha igal juhul
	end subroutine crop_psf
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