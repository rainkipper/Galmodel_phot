module integration_module
	use constants_module
	implicit none
    private
    real(rk),dimension(1:3),parameter:: gauss_x=[ 0.0_rk,-sqrt(3.0_rk/5.0_rk),sqrt(3.0_rk/5.0_rk) ] ! gauss points (do not modify)
    real(rk),dimension(1:3),parameter:: gauss_w=[ 8.0_rk/9.0_rk, 5.0_rk/9.0_rk, 5.0_rk/9.0_rk ] ! gauss weights (do not modify)
	real(rk),parameter:: tiny_value=epsilon(1.0)*1.0E3
	!
    integer,parameter:: n_max_adapt = 10               ! maximum number of adaptive steps (if nn not present)
	integer,parameter:: default_funval_storage_size=1000 ! kui suur on algne fun_val storage size
	integer,parameter:: default_infinity_max_steps=50      ! max number of steps when integrating to infinity
	!
	! defineerime funktsiooni vaartuste tapsuse... seda saab hiljem laiendada kui vaja
	type :: accuracy_type
		real(rk) :: diff_rel ! see on suhteline viga
		real(rk) :: diff_abs=-1.0 ! see on absoluutne viga (default ei kasutata!)
	end type accuracy_type
	!
	!=== define type for storing function values
	type,private:: storage_element_type
        real(rk):: x ! point location
        real(rk):: w ! weight
        real(rk):: val ! vaartus
	end type storage_element_type
    type:: function_value_storage_type
        type(storage_element_type),dimension(:),allocatable:: p ! point along integration line
        integer:: n=0 ! number of points in use
	contains
		procedure,private,pass(this):: check_and_fix_alloc
		procedure,private,pass(this):: initialize_storage
		procedure,public,pass(this):: append_storage
		procedure,private:: copy_fun_storage_value
		generic:: assignment(=) => copy_fun_storage_value
    end type function_value_storage_type
	!
    interface integrate
        module procedure integrate_gauss_adaptive, integrate_gauss_adaptive_noarrfun, &
		integrate_gauss_adaptive_inf, integrate_gauss_adaptive_inf_noarrfun
    end interface
	!
	!========================
	! integrate_gauss_adaptive_inf(func,a,acc, nmax,nmin,errcode,incracc,fvalue,ninf)
	! integrate_gauss_adaptive(func,a,b,acc,   nmax,nmin,errcode,incracc,fvalue,log_base)
	!========================
	!
	public:: accuracy_type, function_value_storage_type
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
    function integrate_gauss_adaptive_inf(func,a,acc,nmax,nmin,errcode,incracc,fvalue,ninf) result(res)
        implicit none
		real(rk),parameter:: multiply=log10(10.0_rk)   ! multiplication parameter - increasing the upper bound
		real(rk),parameter:: min_a=1.0_rk              ! minimum value for lower bound (before starting integratin to infinity)
        real(rk),parameter:: mult_a=1.1_rk             ! multiplication parameter for a if first result is zero (a < min_a)
        real(rk),parameter:: add_a=0.1_rk              ! same as mult_a, but added to a value (a > min_a)
		!
		real(rk) :: res
        real(rk),intent(in):: a     ! integration bounds
        real(rk),intent(in):: acc
        integer,intent(in),optional:: nmax            ! max number of adaptive steps
        integer,intent(in),optional:: nmin         ! min number of adaptive steps
		integer,intent(in),optional:: ninf ! nr of steps to infinity (default vt. ylevalt)
        logical,intent(in),optional:: incracc       ! if true, increase accuracy, do second level down
        character(len=*),intent(in),optional:: errcode
        type(function_value_storage_type),intent(inout),optional :: fvalue ! function value, location and weight arrays
		! parameetrid, mis on lokaalselt kasutusel
		type(function_value_storage_type):: fval
		integer:: nadapt,nadaptmin                  ! max and min nr of adaptive steps
		logical:: increase_acc
		integer,dimension(1:3):: fout      ! temporary variable
		real(rk):: aa,bb,resfirst,res1
		type(accuracy_type):: tacc
		integer:: idum
		integer:: infinity_max_steps
        ! define the interface for integration function
        interface
            function func(x) result(res)
				import:: rk
                implicit none
                real(rk),dimension(:),intent(in) :: x      !coordinate points
                real(rk),dimension(size(x))      :: res    !function values
            end function func
        end interface
		! code starts here
        call fval%initialize_storage() ! initsialiseerime sisemise massiivi funktsiooni vaartustele
        if (present(nmax)) then
            nadapt = nmax
        else
            nadapt = n_max_adapt
        end if
        if (.not.present(nmin)) then
            nadaptmin = 2
        else
            nadaptmin = nmin
        end if
        if (.not.present(incracc)) then
            increase_acc=.false.
        else
            increase_acc=incracc
        end if
		if (present(ninf)) then
			infinity_max_steps = ninf
		else
			infinity_max_steps = default_infinity_max_steps
		end if
		!
		tacc%diff_rel = acc
		!
        ! setting the lower bound according to a value
        ! and calculating first estimates
        if (a>0.0_rk) then ! lower bound is positive
            ! setting the upper value bb - sisuliselt suvaliselt valitud vaartused
            aa=log(a)
            if (a<min_a) then
                bb=log(min_a)
            else
                bb=aa+multiply
            end if
            ! calculating integration result (aa to bb)
			call int_gauss_log(res,func,aa,bb,fout,fval) ! need to be calculated before adaptive one
	        call adaptive_gauss(res,func,aa,bb,res,fout,1,tacc,nadapt,nadaptmin,increase_acc,0.0_rk,fval,log_base=.true.)
			!
            ! if result is zero, decrease upper bound bb
            !  -- probably integration function decrease very (too) rapidly
            if ( abs(res)<= epsilon(0.0) ) then
                if (a<min_a) then
                    bb=log(a*mult_a)
                else
                    bb=aa+add_a
                end if
				call int_gauss_log(res,func,aa,bb,fout,fval) ! need to be calculated before adaptive one
		        call adaptive_gauss(res,func,aa,bb,res,fout,1,tacc,nadapt,nadaptmin+1,.true.,0.0_rk,fval,log_base=.true.)
            end if
            ! if res is still zero, try again... probably tuleb kaugemale integreerida!
            if ( abs(res)<= epsilon(0.0) ) then
                do idum=1,30 ! see on suvaline piir, kuhu ei tohiks kunagi jouda
                    bb=bb+log10(3.0)
					call int_gauss_log(res,func,aa,bb,fout,fval) ! need to be calculated before adaptive one
			        call adaptive_gauss(res,func,aa,bb,res,fout,1,tacc,nadapt,nadaptmin+1,.true.,0.0_rk,fval,log_base=.true.)
                    if (res>0.0_rk) exit ! kui on pos tulemus, siis valju
                end do
            end if
            ! -- if res is still zero, then integration result is also zero
        else ! calculating first estimates, if a is nonpositive
            ! minnakse yha suurtemate bb-de poole, kuna nulli kasvavaid funktsioone ei ole nii palju
            aa=a
            bb=min_a/(10**multiply) ! bb is multiply times smaller than min_a
			call int_gauss_lin(res,func,aa,bb,fout,fval) ! need to be calculated before adaptive one
	        call adaptive_gauss(res,func,aa,bb,res,fout,1,tacc,nadapt,nadaptmin,increase_acc,0.0_rk,fval,log_base=.false.)
            if ( abs(res)<= epsilon(0.0) ) then
                do idum=1,30 ! see on suvaline piir, kuhu ei tohiks kunagi jouda
                    bb=bb*3
					call int_gauss_lin(res,func,aa,bb,fout,fval) ! need to be calculated before adaptive one
			        call adaptive_gauss(res,func,aa,bb,res,fout,1,tacc,nadapt,nadaptmin+1,.true.,0.0_rk,fval,log_base=.false.)
                    if (res>0.0_rk) exit
                end do
            end if
            bb=log(bb) ! see on vaartust, kust hakatakse lopmatusesse integreerima
        end if
		!
		!--------------------------------------------
        ! integration to infinity: alates bb-st
        do idum=1,infinity_max_steps
            resfirst=res ! set the current best value for integration result
            aa=bb
            bb=bb+multiply
            ! calculate integration from aa to bb
			call int_gauss_log(res1,func,aa,bb,fout,fval) ! need to be calculated before adaptive one
	        call adaptive_gauss(res1,func,aa,bb,res,fout,1,tacc,nadapt,nadaptmin,increase_acc,resfirst,fval,log_base=.true.)
            res=res+res1 ! update the current best value
			!
            ! if current best value is accurate enough (according to acc value) exit the cycle
            if (abs(res)*tacc%diff_rel > res1*(10**multiply)) exit
            if ( abs(res1)<= epsilon(0.0) ) exit ! if convergence achieved, exit
        end do
        ! print warning if integration is incomplete
        if (idum==infinity_max_steps) then
            print*, "warning: integrate to infinity incomplete!", idum
            if (present(errcode)) then
                print*, "   ", errcode
            end if
        end if
		!
		! anname valjundisse fvalue kui noutud
        if (present(fvalue)) then
            fvalue = fval
        end if
    end function integrate_gauss_adaptive_inf
	!
    function integrate_gauss_adaptive_inf_noarrfun(func,a,acc,nmax,nmin,errcode,incracc,fvalue,ninf) result(res)
        implicit none
		real(rk) :: res
        real(rk),intent(in):: a     ! integration bounds
        real(rk),intent(in):: acc
        integer,intent(in),optional:: nmax            ! max number of adaptive steps
        integer,intent(in),optional:: nmin         ! min number of adaptive steps
		integer,intent(in),optional:: ninf ! nr of steps to infinity (default vt. ylevalt)
        logical,intent(in),optional:: incracc       ! if true, increase accuracy, do second level down
        character(len=*),intent(in),optional:: errcode
        type(function_value_storage_type),intent(inout),optional :: fvalue ! function value, location and weight arrays
		! parameetrid, mis on lokaalselt kasutusel
		integer:: nadapt,nadaptmin                  ! max and min nr of adaptive steps
		logical:: increase_acc
		integer:: infinity_max_steps
		character(:),allocatable:: ecode
        ! define the interface for integration function
        interface
            function func(x) result(res)
				import:: rk
                implicit none
                real(rk),intent(in) :: x      !coordinate points
                real(rk)      :: res    !function values
            end function func
        end interface
		! code starts here
        if (present(nmax)) then
            nadapt = nmax
        else
            nadapt = n_max_adapt
        end if
        if (.not.present(nmin)) then
            nadaptmin = 2
        else
            nadaptmin = nmin
        end if
        if (.not.present(incracc)) then
            increase_acc=.false.
        else
            increase_acc=incracc
        end if
		if (present(ninf)) then
			infinity_max_steps = ninf
		else
			infinity_max_steps = default_infinity_max_steps
		end if
		if (present(errcode)) then
			ecode=errcode
		else
			ecode = 'kein error message (inf)!'
		end if
		if (present(ninf)) then
			infinity_max_steps = ninf
		else
			infinity_max_steps = default_infinity_max_steps
		end if
		!
		if (present(fvalue)) then
			res = integrate_gauss_adaptive_inf(func_arr,a,acc,nadapt,nadaptmin,ecode,increase_acc,fvalue,infinity_max_steps)
		else
			res = integrate_gauss_adaptive_inf(func_arr,a,acc,nadapt,nadaptmin,ecode,increase_acc,ninf=infinity_max_steps)
		end if
    contains
        Function func_arr(x2) result(res2)
          Implicit None
          real(rk),dimension(:),intent(in) :: x2      !coordinate points
          real(rk),dimension(size(x2))     :: res2    !function values
          integer:: i
          do i=1,size(x2)
            res2(i)=func(x2(i))
          end do
        End Function func_arr
    end function integrate_gauss_adaptive_inf_noarrfun
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
    function integrate_gauss_adaptive(func,a,b,acc,nmax,nmin,errcode,incracc,fvalue,log_base) result(res)
        implicit none
		real(rk) :: res
        real(rk),intent(in):: a,b     ! integration bounds
        real(rk),intent(in):: acc
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
		type(accuracy_type):: tacc
        ! define the interface for integration function
        interface
            function func(x) result(res)
				import:: rk
                implicit none
                real(rk),dimension(:),intent(in) :: x      !coordinate points
                real(rk),dimension(size(x))      :: res    !function values
            end function func
        end interface
		! code starts here
        call fval%initialize_storage() ! initsialiseerime sisemise massiivi funktsiooni vaartustele
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
                print*, "warning: invalid integer bounds (adaptive)!",a,b, blog
            else
                print*, "warning: invalid integer bounds (adaptive)!", errcode, a,b,blog
            end if
            print*, "         : setting result to zero!!!"
            res = 0.0_rk
            return
        else if ( abs(a-b)<= epsilon(0.0) ) then
            res = 0.0_rk
            return
        end if
		!
		tacc%diff_rel = acc
		!
		! start the integrator
		if (blog) then
			aa=log(a); bb=log(b)
			call int_gauss_log(res,func,aa,bb,fout,fval) ! need to be calculated before adaptive one
		else
			aa=a; bb=b
			call int_gauss_lin(res,func,aa,bb,fout,fval) ! need to be calculated before adaptive one
		end if
        call adaptive_gauss(res,func,aa,bb,res,fout,1,tacc,nadapt,nadaptmin,increase_acc,0.0_rk,fval,blog)
		!
		! anname valjundisse fvalue kui noutud
        if (present(fvalue)) then
            fvalue = fval
        end if
    end function integrate_gauss_adaptive
	!
    function integrate_gauss_adaptive_noarrfun(func,a,b,acc,nmax,nmin,errcode,incracc,fvalue,log_base) result(res)
        implicit none
		real(rk) :: res
        real(rk),intent(in):: a,b     ! integration bounds
        real(rk),intent(in):: acc
        integer,intent(in),optional:: nmax            ! max number of adaptive steps
        integer,intent(in),optional:: nmin         ! min number of adaptive steps
        logical,intent(in),optional:: incracc       ! if true, increase accuracy, do second level down
		logical,intent(in),optional:: log_base      ! logaritmic base, default is linear
        character(len=*),intent(in),optional:: errcode
        type(function_value_storage_type),intent(inout),optional :: fvalue ! function value, location and weight arrays
		! parameetrid, mis on lokaalselt kasutusel
		integer:: nadapt,nadaptmin                  ! max and min nr of adaptive steps
		logical:: increase_acc,blog
		character(:),allocatable:: ecode
        ! define the interface for integration function
        interface
            function func(x) result(res)
				import:: rk
                implicit none
                real(rk),intent(in) :: x      !coordinate points
                real(rk)      :: res    !function values
            end function func
        end interface
		! code starts here
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
		if (present(errcode)) then
			ecode=errcode
		else
			ecode = 'kein error message!'
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
        !
		if (present(fvalue)) then
			res = integrate_gauss_adaptive(func_arr,a,b,acc,nadapt,nadaptmin,ecode,increase_acc,fvalue,blog)
		else
			res = integrate_gauss_adaptive(func_arr,a,b,acc,nadapt,nadaptmin,ecode,increase_acc,log_base=blog)
		end if
    contains
        Function func_arr(x2) result(res2)
          Implicit None
          real(rk),dimension(:),intent(in) :: x2      !coordinate points
          real(rk),dimension(size(x2))     :: res2    !function values
          integer:: i
          do i=1,size(x2)
            res2(i)=func(x2(i))
          end do
        End Function func_arr
    end function integrate_gauss_adaptive_noarrfun
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
	! @param resadd - add integral results from previous interval
	!                ii>1 absolute difference
    ! @param resfirst - first estimates for integration - used in accuracy testing
	! @param log_base - use logarithmic base
    ! @result res - integration result
    recursive subroutine adaptive_gauss(res,func,a,b,resold,funin,ii,acc,nn,nnmin,incracc,resadd,fval,log_base)
        implicit none
		! salvestame massiivi kiireks kasutuseks
        real(rk),dimension(1:10),parameter:: iast3=(/1.0d0,3.0d0,9.0d0,27.0d0,81.0d0,243.0d0,&
        729.0d0,2187.0d0,6561.0d0,19683.0d0/)
		real(rk),intent(out):: res
		real(rk),intent(in):: a,b
		real(rk),intent(in):: resold, resadd
		integer,dimension(1:3),intent(in):: funin
		integer,intent(in) :: ii
		class(accuracy_type),intent(in):: acc
		integer,intent(in) :: nn,nnmin
		logical,intent(in):: incracc ! do we increase accuracy -- for sanity check and safety
		type(function_value_storage_type),intent(inout) :: fval ! function value, location and weight arrays
		logical,intent(in):: log_base
		! edasi on ainult lokaalsed muutujad
		real(rk)                   :: t1,t2 ! temporary variables for boundaries, between a..b
        integer,dimension(1:3,1:3):: funout ! sisaldab funktsiooni vaartusi
        real(rk),dimension(1:4)    :: cout ! sisaldab funktsiooni rajasid kolme piirkonna jaoks
        real(rk),dimension(1:3)    :: resc,resca
		real(rk) :: absdiff, iastdum ! absolute difference
        logical:: acctest,doacc ! dotest - teeb teise testi
        ! define the interface for integration function
        interface
            function func(x) result(res)
				import:: rk
                implicit none
                real(rk),dimension(:),intent(in) :: x      !coordinate points
                real(rk),dimension(size(x))      :: res    !function values
            end function func
        end interface
        ! cout on boundaries for subintegrals
        t1=(2.0d0*gauss_x(3)-1.0d0)*( (b-a)*0.5d0 )
        t2=(a+b)*0.5d0
        cout(1)=a; cout(2)=t2-t1; cout(3)=t2+t1; cout(4)=b
        funout(1,:)=funin(:) ! anname edasi keskmise vaartuse igal pool
		call fval%check_and_fix_alloc(6) ! vaatame, kas fval on piisavalt suur massiiv
        ! arvutame funktsiooni vaartused kolmes alamloigus
		if (log_base) then
	        call int_gauss_log_wing(resc(1),func,funin(1),cout(2),cout(3),funout(2:3,1),fval)
	        call int_gauss_log_wing(resc(2),func,funin(2),cout(1),cout(2),funout(2:3,2),fval)
	        call int_gauss_log_wing(resc(3),func,funin(3),cout(3),cout(4),funout(2:3,3),fval)
		else
	        call int_gauss_lin_wing(resc(1),func,funin(1),cout(2),cout(3),funout(2:3,1),fval)
	        call int_gauss_lin_wing(resc(2),func,funin(2),cout(1),cout(2),funout(2:3,2),fval)
	        call int_gauss_lin_wing(resc(3),func,funin(3),cout(3),cout(4),funout(2:3,3),fval)
		end if
		res = sum(resc)
		!
        if (ii==1) then
            absdiff=abs((resadd+res)*acc%diff_rel)
        else
            absdiff=resadd
        end if
        if (ii<=10) then
            iastdum=iast3(ii)
        else
            iastdum=3**(ii-1)
        end if
        !
        ! checking accuracy ja jooksutame uuesti, kui vaja
		doacc=incracc
		if (ii==nn) then ! joudsime maksimum sygavuseni.. ei arvuta edasi
			acctest=.true.
			doacc=.false.
		else
			! koigepealt sanity check
			if (abs(res-resold)<tiny_value) then
				acctest=.true.
			else
				! siin on nyyd tapsuse test
				! ... natuke allpool on samuti tapsuse test
				acctest = abs(res-resold)<abs(res)*acc%diff_rel
				acctest = acctest .or. abs(res)*iastdum<absdiff
			end if
		end if
        if ( doacc .and. acctest ) then
            acctest=.false.
            if (ii>=nnmin) doacc=.false. ! kui ei ole minimaalset sygavust saavutatud, siis ei muuda
        end if
		!
        if ( .not.acctest .or. ii<nnmin ) then
            call adaptive_gauss(resca(1), func,cout(2),cout(3),resc(1),funout(1:3,1),ii+1,acc,nn,nnmin,doacc,absdiff,fval,log_base)
			!
			! jatame aarmised arvutamata, kui tapsus on piisavalt suur
			if (.not.incracc .and. ii>1) then
				res = resca(1) + sum(resc(2:3))
				if ( abs(res-resold)<abs(res)*acc%diff_rel ) then
					return
				end if
			end if
			!
            call adaptive_gauss(resca(2), func,cout(1),cout(2),resc(2),funout(1:3,2),ii+1,acc,nn,nnmin,doacc,absdiff,fval,log_base)
            call adaptive_gauss(resca(3), func,cout(3),cout(4),resc(3),funout(1:3,3),ii+1,acc,nn,nnmin,doacc,absdiff,fval,log_base)
			res = sum(resca)
			!
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
        real(rk),intent(in)                :: a,b
        integer,intent(in):: midp
        real(rk),intent(out):: res
! 		integer:: i
        ! define the interface for integration function
        interface
            function func(x) result(res)
				import:: rk
                implicit none
                real(rk),dimension(:),intent(in) :: x      !coordinate points
                real(rk),dimension(size(x))      :: res    !function values
            end function func
        end interface
        fval%n=fval%n+2
        fval%p(fval%n-1:fval%n)%x = ( (b-a)*(gauss_x(2:3)) + a+b )*0.5_rk
		fval%p(fval%n-1:fval%n)%val = func( fval%p(fval%n-1:fval%n)%x )
		fval%p(fval%n-1:fval%n)%w = ( gauss_w(2:3) ) * ( (b-a)*0.5_rk )
		res=sum( fval%p(fval%n-1:fval%n)%val * fval%p(fval%n-1:fval%n)%w )
		! fix midpoint weight
        fval%p(midp)%w=( gauss_w(1) ) * ( (b-a)*0.5_rk )
		res = res + fval%p(midp)%val * fval%p(midp)%w
		pout(1:2)=[fval%n-1,fval%n]
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
        real(rk),dimension(1:3):: gpoints
        integer,dimension(1:2),intent(out):: pout
        type(function_value_storage_type),intent(inout) :: fval ! function value, location and weight arrays
        real(rk),intent(in)                :: a,b
        integer,intent(in):: midp
        real(rk),intent(out):: res
! 		integer:: i
        interface
            function func(x) result(res)
				import:: rk
                implicit none
                real(rk),dimension(:),intent(in) :: x      !coordinate points
                real(rk),dimension(size(x))      :: res    !function values
            end function func
        end interface
        gpoints(:)=exp( ((b-a)*(gauss_x(:))+a+b)*0.5_rk )
        fval%n=fval%n+2
        fval%p(fval%n-1:fval%n)%x = gpoints(2:3)
		fval%p(fval%n-1:fval%n)%val = func(gpoints(2:3))
		fval%p(fval%n-1:fval%n)%w = ( gpoints(2:3)*gauss_w(2:3) ) * ( (b-a)*0.5_rk )
		res=sum( fval%p(fval%n-1:fval%n)%val * fval%p(fval%n-1:fval%n)%w )
		! fix midpoint weight
        fval%p(midp)%w = ( gpoints(1)*gauss_w(1) ) * ( (b-a)*0.5_rk )
		res = res + fval%p(midp)%val * fval%p(midp)%w
		pout(1:2)=[fval%n-1,fval%n]
    end subroutine int_gauss_log_wing
	!
    ! @param func - integration function
    ! @param a,b - integration bounds
    ! @param pout - gauss points for output (integer indexes to fval array)
	! @param fval - sisaldab koiki punkte, mis on kasutusel
    ! @result res - integration result
    subroutine int_gauss_lin(res,func,a,b,pout,fval)
        implicit none
        integer,dimension(1:3),intent(out):: pout
        type(function_value_storage_type),intent(inout):: fval ! function value, location and weight arrays
        real(rk),intent(in):: a,b
        real(rk),intent(out):: res
! 		integer:: i
        interface
            function func(x) result(res)
				import:: rk
                implicit none
                real(rk),dimension(:),intent(in) :: x      !coordinate points
                real(rk),dimension(size(x))      :: res    !function values
            end function func
        end interface
        fval%n=fval%n+3
        fval%p(fval%n-2:fval%n)%x = ( (b-a)*(gauss_x(:)) + a+b )*0.5_rk
		fval%p(fval%n-2:fval%n)%val = func( fval%p(fval%n-2:fval%n)%x )
        fval%p(fval%n-2:fval%n)%w = ( gauss_w(:) ) * ( (b-a)*0.5_rk )
		res = sum( fval%p(fval%n-2:fval%n)%val * fval%p(fval%n-2:fval%n)%w )
		pout(1:3)=[fval%n-2,fval%n-1,fval%n] ! indeksid punktidele... valjundisse
    end subroutine int_gauss_lin
	!
    ! @param func - integration function
    ! @param a,b - integration bounds
    ! @param pout - gauss points for output (integer indexes to fval array)
	! @param fval - sisaldab koiki punkte, mis on kasutusel
    ! @result res - integration result
    subroutine int_gauss_log(res,func,a,b,pout,fval)
        implicit none
        integer,dimension(1:3),intent(out):: pout ! integer pointer to output massiiv
        type(function_value_storage_type),intent(inout):: fval ! function value, location and weight arrays
        real(rk),intent(in):: a,b
        real(rk),intent(out):: res
! 		integer:: i
        interface
            function func(x) result(res)
				import:: rk
                implicit none
                real(rk),dimension(:),intent(in) :: x      !coordinate points
                real(rk),dimension(size(x))      :: res    !function values
            end function func
        end interface
        fval%n=fval%n+3
        fval%p(fval%n-2:fval%n)%x = exp(((b-a)*(gauss_x(:))+a+b)*0.5_rk)
		fval%p(fval%n-2:fval%n)%val = func( fval%p(fval%n-2:fval%n)%x )
        fval%p(fval%n-2:fval%n)%w = ( fval%p(fval%n-2:fval%n)%x*gauss_w(:) ) * ( (b-a)*0.5_rk )
		res=sum( fval%p(fval%n-2:fval%n)%val * fval%p(fval%n-2:fval%n)%w )
		pout(1:3)=[fval%n-2,fval%n-1,fval%n] ! indeksid punktidele... valjundisse
		!
    end subroutine int_gauss_log
	!
	!-------------------------------------------------------
	!--- defineerime funktsioonid storage tyybile
	!
    subroutine check_and_fix_alloc(this,npluss)
        implicit none
        class(function_value_storage_type),intent(inout):: this
		integer,intent(in):: npluss ! kui palju on juurde vaja
        type(function_value_storage_type):: fdum
! 		integer:: i
		! kontrollime suurust ning suurendame kui vaja
        if (this%n+npluss>size(this%p)) then
            allocate(fdum%p(1:size(this%p)))
			fdum%p(:)%x=this%p(:)%x
            fdum%p(:)%w=this%p(:)%w
			fdum%p(:)%val=this%p(:)%val
			fdum%n=this%n
            deallocate(this%p)
            allocate(this%p(1:this%n*2)) ! suurendame massiivi 2 korda
			this%p(1:this%n)%x=fdum%p(1:this%n)%x
            this%p(1:this%n)%w=fdum%p(1:this%n)%w
			this%p(1:this%n)%val=fdum%p(1:this%n)%val
            deallocate(fdum%p)
        end if
    end subroutine check_and_fix_alloc
	!
    subroutine initialize_storage(this)
        implicit none
        class(function_value_storage_type),intent(inout):: this
! 		integer:: i
		! esialgne allokeerimine
        if (allocated(this%p)) deallocate(this%p)
        allocate(this%p(1:default_funval_storage_size))
        this%n=0
    end subroutine initialize_storage
	!
    subroutine copy_fun_storage_value(fout,fin)
        implicit none
		class(function_value_storage_type),intent(out):: fout
        class(function_value_storage_type),intent(in):: fin
! 		integer:: i
        if (allocated(fout%p)) then
            deallocate(fout%p)
        end if
        allocate(fout%p(1:fin%n))
        fout%n=fin%n
        fout%p(1:fin%n)%x=fin%p(1:fin%n)%x
        fout%p(1:fin%n)%w=fin%p(1:fin%n)%w
		fout%p(1:fin%n)%val=fin%p(1:fin%n)%val
    end subroutine copy_fun_storage_value
	!
    subroutine append_storage(this,fin)
        implicit none
		class(function_value_storage_type),intent(inout):: this
        class(function_value_storage_type),intent(in):: fin
! 		integer:: i
        call this%check_and_fix_alloc(fin%n)
        this%p(this%n+1:this%n+fin%n)%x=fin%p(1:fin%n)%x
        this%p(this%n+1:this%n+fin%n)%w=fin%p(1:fin%n)%w
		this%p(this%n+1:this%n+fin%n)%val=fin%p(1:fin%n)%val
		this%n=this%n + fin%n
    end subroutine append_storage
	!
end module integration_module