module file_operations_module
	use constants_module
contains
	function read_tabel(loetav_file, N) result(res)
		!vajalik spektrite ja filtrite jm andmefailide sisselugemisel
		implicit none
		character(len=default_character_length), intent(in) :: loetav_file
		integer, intent(in) :: N !veergude arv tabelis
		integer :: iunit, ierr
		integer :: ios
		integer :: i
		integer :: ridu
		real(rk), dimension(:,:), allocatable :: res
		character(len=default_character_length) :: rida, formaat
		logical :: olemas
		
		INQUIRE(FILE = trim(loetav_file), exist = olemas)
		if(.not. olemas) then
			print*, "File", trim(loetav_file), "does not exits"
			stop
		end if
		
		formaat = repeat(" ", default_character_length)
		write(unit = formaat, fmt="(A,I,A)", iostat=ios) "(", N, "F)"

		ios  = 0
		iunit = 40
		ridu = 0
		
		!
		! alguses loeb mitu rida on failis, et allokeerida oige suurusega massiivid
		!
		open(file = trim(loetav_file), unit = iunit, action = "read")
		do while(ios .ge. 0)
			ridu = ridu + 1
			rida = repeat(" ", default_character_length) !m2lu nullimine ymberkasutamiseks
			read(fmt = "(A)", iostat=ios, unit=iunit) rida
		end do
! 		print*, "ridu kokku", ridu, trim(loetav_file)
! 		print*, "laius", N, trim(loetav_file)
		close(unit = iunit)
		!
		! massiiivi asjade lugemine
		!
		allocate(res(1:(ridu-1), 1:(N)))
		iunit = 41
		res = 0.0
		i = 0
		ios = 0
		open(file = trim(loetav_file), unit = iunit, action = "read")
		do while(ios .ge. 0)
			i = i+1
			read(fmt = *, iostat=ios, unit=iunit) res(i,:)
		end do
		
! 		stop
		
		close(unit = iunit)
		!
		! ==== populatsiooni muutujate tegemine
		!
! 		do i=1,5
! 			pops(i)%population_name = repeat(" ", default_character_length)
! 			allocate(pops(i)%spec%x(1:size(res, 1)))
! 			allocate(pops(i)%spec%f(1:size(res, 1)))
! 			pops(i)%spec%x = res(:,1)
! 			pops(i)%spec%f = res(:,i+1)
! 		end do
	end function read_tabel
    SUBROUTINE read_fits_to_matrix(filename,pix)
      implicit none
  		integer,parameter  :: rsp = kind(1.0)
      integer,parameter  :: rdp = kind(1.0D0)  !RealKind: double
      INTEGER                                     :: status, unit, blocksize, rw, nfound, group, fpixel
      CHARACTER(len=*),intent(in):: filename ! fits file name
      REAL(rk),DIMENSION(:,:),ALLOCATABLE,intent(out):: pix
      !REAL(rk),DIMENSION(:,:),ALLOCATABLE:: dpix
      INTEGER,dimension(2):: pix_dim
      real(rk)                                    :: nullval
      logical                                     :: anynull
	  logical :: kas_olemas
      status = 0
	  
	  kas_olemas = .false.
	  inquire( file=trim(filename), exist=kas_olemas )
	  if(.not.kas_olemas) then
		  print*, "File does not exist: ", trim(filename)
		  stop
	  end if
      call ftgiou(unit, status) ! get io number for fits file....seda ei pea tegelikult kasutama...voib ka ise maarata

      rw = 0 ! 0-readonly, 1-readwrite
      ! blocksize -- ylearune asi...sellele ei tasu tahelepanu poorata.
      call ftopen(unit, filename, rw, blocksize, status)
      ! Get a sequence of numbered keyword values:
      call ftgknj(unit,'NAXIS', 1, 2, pix_dim, nfound, status)
      ! allocate memory
      ALLOCATE(pix(1:pix_dim(1),1:pix_dim(2)))
      !ALLOCATE( dpix( 0:pix_dim(1)-1,0:pix_dim(2)-1 ) )
      ! read values
      group=1
      fpixel=1
      nullval = -999.0
      IF (rk==rdp) THEN
         CALL ftgpvd(unit, group, fpixel, SIZE(pix), nullval, pix, anynull, status)
      else if (rk==rsp) then
         call ftgpve(unit, group, fpixel, size(pix), nullval, pix, anynull, status)
      else
         PRINT*, "ERROR: no subroutine fitsop: ftgpvx (rk do not mach!)"
      END IF
      !pix=dpix
      ! call ftgkyd(unit, 'delta', delta, comment, status)
      ! call ftgkyd(unit, 'scale', scale, comment, status)
      ! call ftgkyd(unit, 'x0', z_corner(1), comment, status)
      ! call ftgkyd(unit, 'y0', z_corner(2), comment, status)
      ! call ftgkyd(unit, 'z0', z_corner(3), comment, status)
      ! call ftgkyd(unit, 'margin', margin, comment, status)
      ! if (present(rank)) call ftgkyj(unit, 'rank', rank, comment, status)
      ! if (present(in_file)) call ftgkls(unit, 'dat file', in_file, comment, status)
      ! if (present(wt_file)) call ftgkls(unit, 'wts file', wt_file, comment, status)
      call ftclos(unit, status)
      call ftfiou(unit, status)
    END SUBROUTINE read_fits_to_matrix
    SUBROUTINE write_matrix_to_fits(pix,filename)
      use ifport
      implicit none
  	integer,parameter  :: rsp = kind(1.0)
      integer,parameter  :: rdp = kind(1.0D0)  !RealKind: double
      CHARACTER(len=*),intent(in):: filename ! fits file name
      REAL(rk),DIMENSION(:,:),intent(in):: pix
      REAL(rsp),DIMENSION(:,:),allocatable:: pix_s
      INTEGER                                     :: status, unit, blocksize, bitpix, group, fpixel,rw, errnum
      INTEGER,dimension(2):: pix_dim
      logical                                     :: simple, extend
      ALLOCATE(pix_s(1:SIZE(pix,dim=1),1:SIZE(pix,dim=2)))
      pix_s=REAL(pix,kind=rsp)
      status = 0
      call ftgiou(unit, status)
      blocksize = 1
      status = SYSTEM("rm -f "//filename)
      If (status .eq. -1) then
         errnum = ierrno( )
         print *, 'Error ', errnum
      end if
      status=0

      call ftinit(unit, filename, blocksize, status)
      if (status==105) then ! fail eksisteerib
         rw=1
         status=0
         call ftopen(unit, filename, rw, blocksize, status)
         !call ftdopn(unit, filename, rw, status)
      end if
      simple = .true.
      !        bitpix = -64                                           ! double
      bitpix = -32                                                    ! single
      extend = .true.
      pix_dim(1)=SIZE(pix,dim=1)
      pix_dim(2)=SIZE(pix,dim=2)
      call ftphpr(unit, simple, bitpix, size(pix_dim), pix_dim, 0, 1, extend, status)
      group = 1
      fpixel = 1
      IF (rk==rdp) THEN
         !call ftpprd(unit, group, fpixel, size(pix), pix, status)        ! double (Peab vastama pix tyybile!!!)
         call ftppre(unit, group, fpixel, size(pix_s), pix_s, status)        ! single
      else if (rk==rsp) then
         call ftppre(unit, group, fpixel, size(pix), pix, status)        ! single
         PRINT*, "ERROR: no subroutine fitsop: ftpprx (rk do not mach!)"
      END IF

      ! call ftpkyd(unit, 'delta', delta, 5, 'Grid cell size', status)
      ! call ftpkyd(unit, 'scale', scale, 5, 'B3 kernel radius', status)
      ! call ftpkyd(unit, 'margin', margin, 5, 'Margin', status)
      ! call ftpkyd(unit, 'x0', z_corner(1), 10, 'Grid corner location', status)
      ! call ftpkyd(unit, 'y0', z_corner(2), 10, 'Grid corner location', status)
      ! call ftpkyd(unit, 'z0', z_corner(3), 10, 'Grid corner location', status)
      ! if (present(rank)) call ftpkyj(unit, 'rank', rank, 'A trous rank', status)
      ! if (present(in_file)) call ftpkls(unit, 'dat file', in_file, ' ', status)
      ! if (present(wt_file)) call ftpkls(unit, 'wts file', wt_file, ' ', status)

      call ftclos(unit, status)
      call ftfiou(unit, status)
      deallocate(pix_s)
    end subroutine write_matrix_to_fits
end module file_operations_module