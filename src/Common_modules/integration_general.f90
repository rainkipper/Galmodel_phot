module integration_general
	use realkind
	use base_value_type_module
	implicit none
    private
    real(rk),dimension(1:3),parameter:: gauss_x=[ 0.0_rk,-sqrt(3.0_rk/5.0_rk),sqrt(3.0_rk/5.0_rk) ] ! gauss points (do not modify)
    real(rk),dimension(1:3),parameter:: gauss_w=[ 8.0_rk/9.0_rk, 5.0_rk/9.0_rk, 5.0_rk/9.0_rk ] ! gauss weights (do not modify)
	!
    integer,parameter:: n_max_adapt = 10               ! maximum number of adaptive steps (if nn not present)
	integer,parameter:: default_funval_storage_size=1000 ! kui suur on algne fun_val storage size
	!
	!=== define type for storing function values
	type,private:: storage_element_type
        real(rk):: x ! point location
        real(rk):: w ! weight
        class(base_value_type),allocatable:: val ! vaartus
	end type storage_element_type
    type, public:: function_value_storage_type
        type(storage_element_type),dimension(:),allocatable:: p ! point along integration line
        integer:: n=0 ! number of points in use
	contains
		procedure,private,pass(this):: check_and_fix_alloc
		procedure,private,pass(this):: initialize_storage
		procedure,private:: copy_fun_storage_value
		generic:: assignment(=) => copy_fun_storage_value
		final:: finalize_storage
    end type function_value_storage_type
	!
    interface integrate
        module procedure integrate_gauss_adaptive
    end interface
	!
	public:: integrate
	!
contains
	!
    ! @param func - function for integration:  func(x)
    ! @param a - lower bound
	! @param b - upper bound
    ! @param acc - relative accuracy
    ! @param nmax,nmin - maximum and minimum number of adaptive steps
	! @param errcode - error message...hea probleemide lokaliseerimiseks
	! @param incracc - increase accuracy by one level... for testing and sanity check
	! @param fvalue - valjastab vaartused, kus funktsioon on arvutatud, koos kaaludega
	! @param log_base - calculate in logaritmic base, default is linear
    ! @result res - integration result
    !
    subroutine integrate_gauss_adaptive(res,func,a,b,acc,nmax,nmin,errcode,incracc,fvalue,log_base)
        implicit none
		class(base_value_type),intent(inout):: res
        real(rk),intent(in):: a,b     ! integration bounds
        class(accuracy_base_type),intent(in):: acc
        integer,intent(in),optional:: nmax            ! max number of adaptive steps
        integer,intent(in),optional:: nmin         ! min number of adaptive steps
        logical,intent(in),optional:: incracc       ! if true, increase accuracy, do second level down
		logical,intent(in),optional:: log_base      ! logaritmic base, default is linear
        character(len=*),intent(in),optional:: errcode
        type(function_value_storage_type),intent(inout),optional :: fvalue ! function value, location and weight arrays
		! parameetrid, mis on lokaalselt kasutusel
		type(function_value_storage_type):: fval
		integer:: nadapt,nadaptmin                  ! max and min nr of adaptive steps
		logical:: increase_acc,blog
		integer,dimension(1:3):: fout      ! temporary variable
		real(rk):: aa,bb
        ! define the interface for integration function
        interface
            function func(x) result(res)
				import :: rk, base_value_type
                implicit none
                real(rk),intent(in) :: x !coordinate point
                class(base_value_type),allocatable :: res    !function value
            end function func
        end interface
		! code starts here
        call fval%initialize_storage(res) ! initsialiseerime sisemise massiivi funktsiooni vaartustele
        if (present(nmax)) then
            nadapt = nmax
        else
            nadapt = n_max_adapt
        end if
        if (.not.present(nmin)) then
            nadaptmin = 1
        else
            nadaptmin = nmin
        end if
        if (.not.present(incracc)) then
            increase_acc=.false.
        else
            increase_acc=incracc
        end if
        if (.not.present(log_base)) then
            blog=.false.
        else
            blog=log_base
        end if
        ! test the integer bounds
        if (a>b .or. (blog .and. a<=0.0_rk) ) then
            if (.not.present(errcode)) then
                print*, "warning: invalid integer bounds (general)!",a,b, blog
            else
                print*, "warning: invalid integer bounds (general)!", errcode, a,b,blog
            end if
            print*, "         : setting result to zero!!!"
            call res%set_null()
            return
        else if ( abs(a-b)<= epsilon(0.0) ) then
            call res%set_null()
            return
        end if
		! start the integrator
		if (blog) then
			aa=log(a); bb=log(b)
			call int_gauss_log(res,func,aa,bb,fout,fval) ! need to be calculated before adaptive one
		else
			aa=a; bb=b
			call int_gauss_lin(res,func,aa,bb,fout,fval) ! need to be calculated before adaptive one
		end if
        call adaptive_gauss(res,func,aa,bb,res,fout,1,acc,nadapt,nadaptmin,increase_acc,fval,present(fvalue),blog)
		!
		! anname valjundisse fvalue kui noutud
        if (present(fvalue)) then
            fvalue = fval
        end if
    end subroutine integrate_gauss_adaptive
	!
    ! @param func - integration function: fun(x)
    ! @param a,b - integration bounds
    ! @param resold - integration result in previous step: this value is used for testing accuracy
    ! @param funin - gauss points from previous step: these points are middle gauss points
    !                in this step, calculate only two gauss points (gauss wings)
    ! @param ii - index of adaptive step
    ! @param acc - relative accuracy
    ! @param nn,nnmin - max and min nr of adaptive steps
	! @param incracc - increase accuracy one level for testing
    ! @param resfirst - first estimates for integration - used in accuracy testing
	! @param store_fval - store fval vaartus for output?
	! @param log_base - use logarithmic base
    ! @result res - integration result
    recursive subroutine adaptive_gauss(res,func,a,b,resold,funin,ii,acc,nn,nnmin,incracc,fval,store_fval,log_base)
        implicit none
		class(base_value_type),intent(inout):: res
		real(rk),intent(in):: a,b
		class(base_value_type),intent(in):: resold
		integer,dimension(1:3),intent(in):: funin
		integer,intent(in) :: ii
		class(accuracy_base_type),intent(in):: acc
		integer,intent(in) :: nn,nnmin
		logical,intent(in):: incracc ! do we increase accuracy -- for sanity check and safety
		type(function_value_storage_type),intent(inout) :: fval ! function value, location and weight arrays
		logical,intent(in):: store_fval, log_base
		! edasi on ainult lokaalsed muutujad
		real(rk)                   :: t1,t2 ! temporary variables for boundaries, between a..b
        integer,dimension(1:3,1:3):: funout ! sisaldab funktsiooni vaartusi
        real(rk),dimension(1:4)    :: cout ! sisaldab funktsiooni rajasid kolme piirkonna jaoks
        class(base_value_type),allocatable    :: resc1,resc2,resc3 !
		class(base_value_type),allocatable    :: resc1a,resc2a,resc3a
        logical:: acctest,doacc ! dotest - teeb teise testi
		integer:: i,k
        ! define the interface for integration function
        interface
            function func(x) result(res)
				import :: rk, base_value_type
                implicit none
                real(rk),intent(in) :: x !coordinate point
                class(base_value_type),allocatable :: res    !function value
            end function func
        end interface
        ! cout on boundaries for subintegrals
        t1=(2.0d0*gauss_x(3)-1.0d0)*( (b-a)*0.5d0 )
        t2=(a+b)*0.5d0
        cout(1)=a; cout(2)=t2-t1; cout(3)=t2+t1; cout(4)=b
        funout(1,:)=funin(:) ! anname edasi keskmise vaartuse igal pool
		call fval%check_and_fix_alloc(res) ! vaatame, kas fval on piisavalt suur massiiv
        ! ylejaanud punktid tulevad jargmistest ridadest
        ! arvutame funktsiooni vaartused kolmes alamloigus
		call res%alloc_same_type(resc1)
		call res%alloc_same_type(resc2)
		call res%alloc_same_type(resc3)
		if (log_base) then
	        call int_gauss_log_wing(resc1,func,funin(1),cout(2),cout(3),funout(2:3,1),fval)
	        call int_gauss_log_wing(resc2,func,funin(2),cout(1),cout(2),funout(2:3,2),fval)
	        call int_gauss_log_wing(resc3,func,funin(3),cout(3),cout(4),funout(2:3,3),fval)
		else
	        call int_gauss_lin_wing(resc1,func,funin(1),cout(2),cout(3),funout(2:3,1),fval)
	        call int_gauss_lin_wing(resc2,func,funin(2),cout(1),cout(2),funout(2:3,2),fval)
	        call int_gauss_lin_wing(resc3,func,funin(3),cout(3),cout(4),funout(2:3,3),fval)
		end if
		call res%set_null()
		call res%add(resc1)
		call res%add(resc2)
		call res%add(resc3)
        !
        ! checking accuracy ja jooksutame uuesti, kui vaja
		doacc=incracc
		if (ii==nn) then ! joudsime maksimum sygavuseni.. ei arvuta edasi
			acctest=.true.
			doacc=.false.
		else
			acctest=res%is_similar_with(resold,acc)
		end if
        if ( doacc .and. acctest ) then
            acctest=.false.
            if (ii>=nnmin) doacc=.false. ! kui ei ole minimaalset sygavust saavutatud, siis ei muuda
        end if
		!
        if ( .not.acctest .or. ii<nnmin ) then
			call res%alloc_same_type(resc1a)
			call res%alloc_same_type(resc2a)
			call res%alloc_same_type(resc3a)
            call adaptive_gauss(resc1a, func,cout(2),cout(3),resc1,funout(1:3,1),ii+1,acc,nn,nnmin,doacc,fval,store_fval,log_base)
			!
			! jatame aarmised arvutamata, kui tapsus on piisavalt suur
			if (.not.incracc .and. ii>1) then
				call res%set_null()
				call res%add(resc1a)
				call res%add(resc2)
				call res%add(resc3)
				if (res%is_similar_with(resold,acc)) then
					return
				end if
			end if
			!
            call adaptive_gauss(resc2a, func,cout(1),cout(2),resc2,funout(1:3,2),ii+1,acc,nn,nnmin,doacc,fval,store_fval,log_base)
            call adaptive_gauss(resc3a, func,cout(3),cout(4),resc3,funout(1:3,3),ii+1,acc,nn,nnmin,doacc,fval,store_fval,log_base)
			call res%set_null()
			call res%add(resc1a)
			call res%add(resc2a)
			call res%add(resc3a)
			!
		else if (.not.store_fval) then
			do k=1,3
				do i=1,3
					call fval%p(funout(i,k))%val%dealloc()
					deallocate(fval%p(funout(i,k))%val)
				end do
			end do
        end if
    end subroutine adaptive_gauss
	!
    ! @param func - integration function
    ! @param midp - middle point value
    ! @param a,b - integration bounds
    ! @param pout - gauss points (wings) for output (integer indexes to fval array)
	! @param fval - sisaldab koiki punkte, mis on kasutusel
    ! @result res - integration result
    subroutine int_gauss_lin_wing(res,func,midp,a,b,pout,fval)
        implicit none
        integer,dimension(1:2),intent(out):: pout
        type(function_value_storage_type),intent(inout) :: fval ! function value, location and weight arrays
        real(rk),dimension(1:2)            :: gpoints ! temporary variable
        real(rk),intent(in)                :: a,b
        integer,intent(in):: midp
        class(base_value_type),intent(inout):: res
		integer:: i
        ! define the interface for integration function
        interface
            function func(x) result(res)
				import :: rk, base_value_type
                implicit none
                real(rk),intent(in) :: x !coordinate point
                class(base_value_type),allocatable :: res    !function value
            end function func
        end interface
        gpoints(:)=( (b-a)*(gauss_x(2:3)) + a+b )*0.5_rk
        fval%n=fval%n+2
        fval%p(fval%n-1:fval%n)%x=gpoints(1:2)
		call fval%p(fval%n-1)%val%copy(func(gpoints(1)))
		call fval%p(fval%n)%val%copy(func(gpoints(2)))
        pout(1:2)=[fval%n-1,fval%n]
        fval%p(fval%n-1:fval%n)%w=( gauss_w(2:3) ) * ( (b-a)*0.5_rk )
		call res%set_null()
		do i=fval%n-1,fval%n
			call res%add( fval%p(i)%val , fval%p(i)%w )
		end do
		! fix midpoint weight
        fval%p(midp)%w=( gauss_w(1) ) * ( (b-a)*0.5_rk )
		call res%add( fval%p(midp)%val , fval%p(midp)%w )
    end subroutine int_gauss_lin_wing
    !
    ! calculating the gauss wings - middle points are input points
    ! using logarithmic scale
    !
    ! @param func - integration function
    ! @param midp - middle point value
    ! @param a,b - integration bounds
    ! @param pout - gauss points (wings) for output (integer indexes to fval array)
	! @param fval - sisaldab koiki punkte, mis on kasutusel
    ! @result res - integration result
    subroutine int_gauss_log_wing(res,func,midp,a,b,pout,fval)
        implicit none
        real(rk),dimension(1:3)            :: gpoints
        integer,dimension(1:2),intent(out):: pout
        type(function_value_storage_type),intent(inout) :: fval ! function value, location and weight arrays
        real(rk),intent(in)                :: a,b
        integer,intent(in):: midp
        class(base_value_type),intent(inout):: res
		integer:: i
        interface
            function func(x) result(res)
				import :: rk, base_value_type
                implicit none
                real(rk),intent(in) :: x !coordinate point
                class(base_value_type),allocatable :: res    !function value
            end function func
        end interface
        gpoints(:)=exp( ((b-a)*(gauss_x(:))+a+b)*0.5_rk )
        fval%n=fval%n+2
        fval%p(fval%n-1:fval%n)%x=gpoints(2:3)
		call fval%p(fval%n-1)%val%copy(func(gpoints(2)))
		call fval%p(fval%n)%val%copy(func(gpoints(3)))
		pout(1:2)=[fval%n-1,fval%n]
        fval%p(fval%n-1:fval%n)%w=( gpoints(2:3)*gauss_w(2:3) ) * ( (b-a)*0.5_rk )
		call res%set_null()
		do i=fval%n-1,fval%n
			call res%add( fval%p(i)%val , fval%p(i)%w )
		end do
		! fix midpoint weight
        fval%p(midp)%w=( gpoints(1)*gauss_w(1) ) * ( (b-a)*0.5_rk )
		call res%add( fval%p(midp)%val , fval%p(midp)%w )
    end subroutine int_gauss_log_wing
	!
    ! @param func - integration function
    ! @param a,b - integration bounds
    ! @param pout - gauss points for output (integer indexes to fval array)
	! @param fval - sisaldab koiki punkte, mis on kasutusel
    ! @result res - integration result
    subroutine int_gauss_lin(res,func,a,b,pout,fval)
        implicit none
        real(rk),dimension(1:3):: gpoints
        integer,dimension(1:3),intent(out):: pout
        type(function_value_storage_type),intent(inout):: fval ! function value, location and weight arrays
        real(rk),intent(in):: a,b
        class(base_value_type),intent(inout):: res
		integer:: i
        interface
            function func(x) result(res)
				import :: rk, base_value_type
                implicit none
                real(rk),intent(in) :: x !coordinate point
                class(base_value_type),allocatable :: res    !function value
            end function func
        end interface
		gpoints(:)=( (b-a)*(gauss_x(:)) + a+b )*0.5_rk
        fval%n=fval%n+3
        fval%p(fval%n-2:fval%n)%x=gpoints(1:3)
		do i=fval%n-2,fval%n
			call fval%p(i)%val%copy( func(fval%p(i)%x) )
		end do
		pout(1:3)=[fval%n-2,fval%n-1,fval%n] ! indeksid punktidele... valjundisse
        fval%p(fval%n-2:fval%n)%w=( gauss_w(:) ) * ( (b-a)*0.5_rk )
		!
		call res%set_null()
		do i=fval%n-2,fval%n
			call res%add( fval%p(i)%val , fval%p(i)%w )
		end do
    end subroutine int_gauss_lin
	!
    ! @param func - integration function
    ! @param a,b - integration bounds
    ! @param pout - gauss points for output (integer indexes to fval array)
	! @param fval - sisaldab koiki punkte, mis on kasutusel
    ! @result res - integration result
    subroutine int_gauss_log(res,func,a,b,pout,fval)
        implicit none
        real(rk),dimension(1:3):: gpoints
        integer,dimension(1:3),intent(out):: pout ! integer pointer to output massiiv
        type(function_value_storage_type),intent(inout):: fval ! function value, location and weight arrays
        real(rk),intent(in):: a,b
        class(base_value_type),intent(inout):: res
		integer:: i
        interface
            function func(x) result(res)
				import :: rk, base_value_type
                implicit none
                real(rk),intent(in) :: x !coordinate point
                class(base_value_type),allocatable :: res    !function value
            end function func
        end interface
        gpoints(:)=exp(((b-a)*(gauss_x(:))+a+b)*0.5_rk)
        fval%n=fval%n+3
        fval%p(fval%n-2:fval%n)%x=gpoints(1:3)
		do i=fval%n-2,fval%n
			call fval%p(i)%val%copy( func(fval%p(i)%x) )
		end do
		pout(1:3)=[fval%n-2,fval%n-1,fval%n] ! indeksid punktidele... valjundisse
        fval%p(fval%n-2:fval%n)%w=( gpoints(:)*gauss_w(:) ) * ( (b-a)*0.5_rk )
		!
		call res%set_null()
		do i=fval%n-2,fval%n
			call res%add( fval%p(i)%val , fval%p(i)%w )
		end do
    end subroutine int_gauss_log
	!
	!-------------------------------------------------------
	!--- defineerime funktsioonid storage tyybile
	!
    subroutine check_and_fix_alloc(this,rtype)
        implicit none
        class(function_value_storage_type),intent(inout):: this
        type(function_value_storage_type):: fdum
		class(base_value_type),intent(in):: rtype ! reference type... seda tyypi allokeerime
		integer:: i
		! kontrollime suurust ning suurendame kui vaja
        if (this%n+6>size(this%p)) then
            allocate(fdum%p(1:size(this%p)))
			fdum%p(:)%x=this%p(:)%x
            fdum%p(:)%w=this%p(:)%w
			fdum%n=this%n
			do i=1,this%n
				if (allocated(this%p(i)%val)) then
					call rtype%alloc_same_type(fdum%p(i)%val,construct=.false.)
					call fdum%p(i)%val%copy( this%p(i)%val )
				end if
			end do
            deallocate(this%p)
            allocate(this%p(1:this%n*2)) ! suurendame massiivi 2 korda
			this%p(1:this%n)%x=fdum%p(1:this%n)%x
            this%p(1:this%n)%w=fdum%p(1:this%n)%w
			do i=1,size(this%p)
				if (i<=this%n .and. .not.allocated(fdum%p(i)%val)) cycle
				call rtype%alloc_same_type(this%p(i)%val,construct=.false.)
				if (i>this%n) cycle
				call this%p(i)%val%copy( fdum%p(i)%val )
			end do
            deallocate(fdum%p)
        end if
    end subroutine check_and_fix_alloc
	!
    subroutine initialize_storage(this,rtype)
        implicit none
        class(function_value_storage_type),intent(inout):: this
		class(base_value_type),intent(in):: rtype ! reference type... seda tyypi allokeerime
		integer:: i
		! esialgne allokeerimine
        if (allocated(this%p)) deallocate(this%p)
        allocate(this%p(1:default_funval_storage_size))
        this%n=0
		do i=1,size(this%p)
			call rtype%alloc_same_type(this%p(i)%val,construct=.false.)
		end do
    end subroutine initialize_storage
	!
    subroutine copy_fun_storage_value(fout,fin)
        implicit none
		class(function_value_storage_type),intent(out):: fout
        class(function_value_storage_type),intent(in):: fin
		integer:: i
        if (allocated(fout%p)) then
            deallocate(fout%p)
        end if
        allocate(fout%p(1:fin%n))
        fout%n=fin%n
        fout%p(1:fin%n)%x=fin%p(1:fin%n)%x
        fout%p(1:fin%n)%w=fin%p(1:fin%n)%w
		do i=1,fin%n
			if (allocated(fin%p(i)%val)) then
				call fin%p(i)%val%alloc_same_type(fout%p(i)%val)
				call fout%p(i)%val%copy(fin%p(i)%val)
			end if
		end do
    end subroutine copy_fun_storage_value
	!
	subroutine finalize_storage(this)
		implicit none
		type(function_value_storage_type):: this
		integer:: i
		if (allocated(this%p)) deallocate(this%p)
	end subroutine finalize_storage
	!
end module integration_general