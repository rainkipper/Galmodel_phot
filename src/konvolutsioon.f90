module konvolutsioon_module
    use constants_module
!     use images
!     implicit none
!     real(rk), dimension(:,:), allocatable :: conv
    
!     call convolve("mudelpilt.fits","esimene_mudelpilt.fits",conv)
!     call write_matrix_to_fits(conv, "a.fits")  

contains
    ! kolm fourier transform algoritmi -

    ! DFT(array,size,direction) (dft, standard) | 
    ! four1(array,size,dir) (fft, numerical recipes in C) | eeldab size = 2^x
    ! fft_radix_2(array,size,dir) (fft, Cooley-Turkey algorithm) | eeldab ^
    !    
    ! andmete massiivi (1D) struktuuri naide
    ! kaks kompleksarvu Re1+Im1, Re2+Im2
    ! array = x, size = 2    

    ! NB! Tahele panna, et kuigi massiivi tegelik suurus N on 4, siis
    ! meie kasitleme seda kui N/2 ehk 2, kuna kompleksarv koosneb 
    ! reaalosast ja imaginaarosast, mis asetsevad korvuti

    ! Re1 Im1 Re2 Im2
    ! i=1   2   3   4
    ! 
    ! read_fits_to_matrix(filename,array) 
    ! makepadded(a1_in, a1_in_size, size_after_padding, a1_out)
    ! add_complex(a1_size, a1_inout)
    ! d2_fft(a1_inout, a1_size, direction) - 2D fft
    
    ! main convolution exec flow
    subroutine convolve(pix_img, pix_kern, convolved)
        implicit none
        real(rk), dimension(:,:), allocatable, intent(in) :: pix_img, pix_kern
        real(rk), intent(inout), dimension(:,:), allocatable :: convolved ! convolved 2d array
        real(rk), dimension(:,:), allocatable :: img_padded, kern_padded
        integer :: n1, n2, dn
        integer :: i, m, pow
!         real(rk), allocatable :: pix_sum(:), out_sum(:)

        ! assume len xdim = len ydim
        n1 = int(sqrt(float(size(pix_img))))
        n2 = int(sqrt(float(size(pix_kern))))

        ! find next upper power of two for padding
        if(n2<n1) then
            dn = n1
        else if(n1<=n2) then
            dn = n2
        end if
        m = 0
        pow = 1
        do while(pow < dn) 
            pow = pow*2
            m = m+1
        end do
        
        allocate(img_padded(1:pow,1:pow))
        allocate(kern_padded(1:pow,1:pow))
            
        !pad with zeroes up to next power of two
        call makepadded(pix_img, n1, pow, img_padded)  
        call makepadded(pix_kern, n2, pow, kern_padded)  
! 		print*, "CP1"
! 		print*, size(img_padded)
! 		print*, size(kern_padded)
! 		print*, "CP2"
        
        ! add complex counterparts to real pixel values
        ! sequence: i (odd) - real
        !           i+1 (even) - img
        call add_complex(pow, img_padded) 
        call add_complex(pow, kern_padded)
               
        ! forward transform
        call d2_fft(img_padded, pow, 1)
        call d2_fft(kern_padded, pow, 1) 
        
       ! multiply in frequency domain
        print *,"! multiply in frequency domain"
        allocate(convolved(1:pow,1:2*pow))
       
        do i=1,pow
            call cmplx_multiply(img_padded(i,:), kern_padded(i,:), convolved(i,:))     
        end do
        
        ! inverse fft
        print *,"! inverse fft"
        call d2_fft(convolved, pow, -1)

        print *,"! fix "
        ! transform image pieces to place (because of frequency shift)        
        call fix_img(convolved, pow)

        print *,"! remove complex"
         ! remove complex counterparts
        call remove_complex(pow, convolved)
        
        print *,"! remove padding"
        ! remove added zeroes
        call remove_padding(convolved, pow, dn)
        
    end subroutine convolve

    subroutine fix_img(img, s)
        implicit none
        real(rk), intent(inout) :: img(:,:)
        integer, intent(in) :: s ! size of 1dim
        integer :: i, j
        
        do i=1,s/2
            do j=1,2*s
                call swap(img(i,j),img(s/2+i,j))
            end do
        end do
       
        do i=1,s
            do j=1,s
                call swap(img(i,j),img(i,s+j))
            end do
        end do
    
        do i=1,s
            do j=1,s
                call swap(img(i,j),img(i,2*s-j))
            end do
        end do
    end subroutine fix_img
    
    subroutine makepadded(pix, original_size, padded_size, pix_padded)
        ! pix - 2d image array
        implicit none
        real(rk), dimension(:,:), intent(in) :: pix
        integer, intent(in) :: original_size, padded_size
        real(rk), dimension(:,:), intent(out) :: pix_padded
!         real(rk) :: temp
        integer :: i, j, padx!, c, ks2

        padx = (padded_size-original_size)/2
        
        do i=1, padded_size
            do j=1, padded_size
                if(j > padded_size-padx .or. j <= padx .or. i <= padx .or. i > padded_size-padx) then
                    pix_padded(i,j) = 0._rk
                else
                    pix_padded(i,j) = pix(i-padx,j-padx)
                end if
            end do
        end do

    end subroutine makepadded
    
    subroutine remove_padding(img, s, r)
        implicit none
        real(rk), intent(inout), allocatable :: img(:,:)
        real(rk), allocatable :: temp(:,:)
        integer, intent(in) :: s, r ! m,n dim size; r original size
        integer :: i,j, dif
    
        dif = (s-r)/2
        
        allocate(temp(1:r,1:r))
        do i=1, r
            do j=1, r
                temp(i,j) = img(dif+i,dif+j)
            end do
        end do

        deallocate(img)
        allocate(img(1:r,1:r))
        do i=1, r
            do j=1, r
                img(i,j) = temp(i,j)
            end do
        end do
       
    end subroutine remove_padding

    subroutine add_complex(n, arr)
        implicit none
        integer, intent(in) :: n
        real(rk), intent(inout), dimension(:,:), allocatable :: arr
        real(rk), dimension(:,:), allocatable :: temp_arr
        integer :: i,j, k !,c

        allocate(temp_arr(1:n,1:2*n))

        do i=1,n
            k = 1
            do j=1,2*n,2            
                temp_arr(i,j) = arr(i,k)
                temp_arr(i,j+1) = 0.0_rk
                k = k+1
            end do
        end do  

        deallocate(arr)
        allocate(arr(1:n,1:2*n))
        do i=1,n
            do j=1,2*n
                arr(i,j) = temp_arr(i,j)
            end do
        end do
       
    end subroutine add_complex

    subroutine remove_complex(n, arr)
        implicit none
        integer, intent(in) :: n
        real(rk), intent(inout), dimension(:,:), allocatable :: arr
        real(rk), dimension(:,:), allocatable :: temp_arr
        integer :: i,j,c 
		integer :: ALLOC_ERR
! 		character(len=:), allocatable :: err
        allocate(temp_arr(n,n))
		print*, "CP1", size(arr, 1), size(arr,2)
        do i=1, n
            c = 2*n-1
            do j=1, 2*n
                temp_arr(i,j) = arr(i,c)
                c = c - 2
            end do
        end do
		print*, "CP2", allocated(arr), allocated(temp_arr)
! 		deallocate(arr, 		 ERRMSG=err)
! 		print*, err
		deallocate(arr, STAT = ALLOC_ERR)
! 		deallocate(arr)
		print*, "CP3", ALLOC_ERR
		stop "koik"
        allocate(arr(1:n,1:n))
        do i=1,n
            do j=1,n
                arr(i,j) = temp_arr(i,j)
            end do
        end do
        print*, "CP4"
    end subroutine remove_complex

    subroutine cmplx_multiply(x, y, out)
        implicit none
        ! len(x) = len(y) = len(out)
        !          
        real(rk), dimension(:), intent(in) :: x,y
        real(rk), dimension(:), intent(inout) :: out
        integer :: i
        
        do i=1,size(x),2
            out(i) = x(i)*y(i)-x(i+1)*y(i+1)
            out(i+1) = x(i)*y(i+1)+x(i+1)*y(i)
        end do

    end subroutine cmplx_multiply

    subroutine swap(a,b)
        implicit none
        real(rk), intent(inout) :: a,b
        real(rk) :: temp
        temp = a
        a = b
        b = temp
        
    end subroutine swap

    subroutine four1(dat, nn, isign) ! data, length of data ( img + real parts)
        implicit none
        ! isign - forward/inverse
        real(rk), dimension(:), intent(inout) :: dat
        integer, intent(in) :: nn, isign
        integer :: n, mmax, m, j, istep, i
        real(rk) :: wtemp, wr, wpr, wpi, wi, theta, tempr, tempi!, z
        
        n = ishft(nn, 1)
        j=1
        do i=1,n,2
            if(j>i) then
                call swap(dat(i),dat(j))
                call swap(dat(j+1),dat(i+1))
            end if
            m = ishft(n,-1)
            do while(m >= 2 .and. j > m)
                j = j-m
                m = ishft(m, -1)
            end do
            j = j+m
        end do
  
        mmax = 2
        do while(n > mmax)
            istep = 2*mmax
            theta = 2*PI/(-isign*mmax)
            wtemp = sin(0.5*theta)
            wpr = -2.0*wtemp*wtemp
            wpi = sin(theta)
            wr = 1.0
            wi = 0.0
            do m=1,mmax,2
                do i=m,n,istep
                    j=i+mmax
                    tempr = wr*dat(j)-wi*dat(j+1)
                    tempi=wr*dat(j+1)+wi*dat(j)
                    dat(j) = dat(i) - tempr
                    dat(j+1) = dat(i+1)-tempi
                    dat(i) = dat(i) + tempr
                    dat(i+1) = dat(i+1) + tempi
                end do
                wtemp = wr    
                wr = wr*wpr-wi*wpi+wr
                wi = wi*wpr+wtemp*wpi+wi
            end do
            mmax=istep
        end do 

        if(isign .eq. -1) then  
            dat = dat/nn
        end if  
        
    end subroutine four1
    
    subroutine d2_fft(arr, m, dir)
        ! 2d dft
        implicit none
        integer, intent(in) :: dir, m ! forward(1)/inverse(-1) transform 
        real(rk), dimension(:,:), intent(inout) :: arr
        real(rk), dimension(:), allocatable :: t
        integer :: i,j,k
        
        !1d ft on row
        do i=1,m
            call fft_radix_2(arr(i,:),m, dir)
        end do
        
        !1d ft on col
        allocate(t(1:2*m))
        
        do i=1,2*m,2
            ! teisendame tulbad i ja i+1 reaks 
            ! t(1) = i(1), t(2) = i+1(1) t(3) = i(2) ...
            k = 1 
            do j=1,2*m, 2
                t(j) = arr(k,i)
                t(j+1) = arr(k,i+1)
                k = k+1
            end do
            
            call fft_radix_2(t, m, dir)

            ! teisendame rea tagasi tulbaks
            k = 1            
            do j=1,2*m,2
                arr(k,i) = t(j)
                arr(k,i+1) = t(j+1)
                k = k+1
            end do

        end do

    end subroutine d2_fft

    subroutine DFT(arr, m, dir)
        implicit none
        integer, intent(in) :: m, dir ! complex numbers size
        real(rk), dimension(:), intent(inout) :: arr
        integer :: i,n,k
        complex, dimension(:), allocatable :: temp, arr_c
!         real(rk) :: cosarg, sinarg

        allocate(arr_c(1:m))
        k = 1
        do i=1,m*2,2
            arr_c(k) = cmplx(arr(i),arr(i+1))
            k = k+1
        end do

        allocate(temp(1:m))  
        do i=1,m
            temp(i) = 0.0
        end do
              
        do n=0,m-1
            do k=0,m-1
                temp(n+1) = temp(n+1) + arr_c(k+1)*exp(cmplx(0,-2*dir*PI/m*k*n))
            end do
        end do
        
        k = 1
        do i=1,2*m,2
            arr(i) = real(temp(k))
            arr(i+1) = aimag(temp(k))
            k = k+1
        end do

        if(dir .eq. -1) then
            arr = arr/m
        end if
        
        deallocate(temp)
    end subroutine DFT
    
    subroutine fft_radix_2(arr,N, dir)
        implicit none
        integer, intent(in) :: N, dir
        real(rk), intent(inout), dimension(:) :: arr   
        integer :: i,k,nn!,m
        complex, dimension(:), allocatable :: Xk, arrToComplex
        
        nn = 2*N
        
        allocate(Xk(1:N))
        allocate(arrToComplex(1:N))

        ! convert from real-img pairs to complex
        k = 1
        do i=1,nn,2
            arrToComplex(k) = CMPLX(arr(i),arr(i+1))
            k = k+1
        end do
        
        Xk = radix2fft(arrToComplex,N, 1, dir)                        
    
        if(dir == -1) then
             Xk = Xk/(N)
        end if

        ! convert from complex to real-img pairs
        k = 1
        do i=1,nn,2
            arr(i) = REAL(Xk(k))
            arr(i+1) = AIMAG(Xk(k))
            k = k+1 
        end do

    end subroutine fft_radix_2

    ! radix-2 dit cooley-turkey algorithm
    recursive function radix2fft(arr, N, s, dir) result (Xk)
        implicit none
        integer :: N, nn, s, dir, i,k !,m
        complex, dimension(:), allocatable :: arr, Xk, t1, t2
        real(rk) :: arg

        allocate(Xk(1:N)) 
        if(N /= 1) then

            nn = N/2
            ! jagame kaheks
            allocate(t1(1:nn))
            allocate(t2(1:nn))
            k = 1
            do i=1,N,2
                t1(k) = arr(i)
                t2(k) = arr(i+1)
                k = k+1
            end do
          
            ! molema poole ft
            t1 = radix2fft(t1,nn,2*s,dir)
            t2 = radix2fft(t2,nn,2*s,dir)
    
            arg = -2*PI/N*dir            
            do k=0, nn-1
                Xk(k+1) = t1(k+1) + exp(cmplx(0,arg*k))*t2(k+1)
                Xk(k+nn+1) = t1(k+1) - exp(cmplx(0,arg*k))*t2(k+1)
            end do

        else
            Xk(1) = arr(1)
        end if 

    end function radix2fft
end module konvolutsioon_module
