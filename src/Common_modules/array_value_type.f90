module array_value_type_module
	use realkind
	use base_value_type_module
	implicit none
	private
	public:: array_value_type, base_value_type, accuracy_base_type
	!
	! kasutab lbound ja ubound, et array algust ja loppu testida
	type, extends(base_value_type) :: array_value_type
		real(rk),dimension(:),allocatable:: v ! vector/array 
	contains
		procedure,pass(this),public:: alloc_same_type => alloc_same_type_tarray
		procedure,pass(this),public:: set_null => set_null_tarray
		procedure,pass(this),public:: dealloc => dealloc_tarray
		procedure,pass(this),public:: is_similar_with => is_similar_with_tarray
		! moned private protseduurid, mis on generic all defineeritud
		procedure,pass(this),private:: construct_from_value => construct_from_value_tarray
		procedure,pass(this),private:: construct_from_params => construct_from_params_tarray
		! protseduurid, mis on vajalikud generic funktsioonide jaoks
		procedure,pass(this),private:: copy_value => copy_value_tarray
		procedure,pass(this),private:: copy_real => copy_real_tarray
		procedure,pass(this),private:: sum_value_and_value => sum_value_and_value_tarray
		procedure,pass(this),private:: sum_value_and_real => sum_value_and_real_tarray
		procedure,pass(this),private:: add_value => add_value_tarray
		procedure,pass(this),private:: add_real => add_real_tarray
		procedure,pass(this),private:: add_value_times_constant => add_value_times_constant_tarray
		procedure,pass(this),private:: multiply_with_constant => multiply_with_constant_tarray
		! extended tyypide ja allokeerivate muutujate korral
		! ... tuleb defineerida final funktsioon
		final:: clean_object
	end type array_value_type
	!
contains
	!
	! define kahe vaartuse tapsuse hindamine
	! isel -- sellega saab erinevaid skeeme implementeerida ning nende vahel valida
	function is_similar_with_tarray(this,ref,acc,isel) result(res)
		implicit none
		class(array_value_type),intent(in):: this
		class(base_value_type),intent(in):: ref
		class(accuracy_base_type),intent(in):: acc
		integer,intent(in),optional:: isel ! sellega saab valida, millist eeskirja kasutatakse!
		logical:: res
		res = this%base_value_type%is_similar_with(ref,acc)
	end function is_similar_with_tarray
	!
	! set null value... paneb ka null parameetri nulliks!
	subroutine set_null_tarray(this)
		implicit none
		class(array_value_type),intent(inout):: this
		call this%base_value_type%set_null()
		if (allocated(this%v)) this%v=0.0_rk
	end subroutine set_null_tarray
	! see on ainult kohahoidja
	subroutine dealloc_tarray(this)
		implicit none
		class(array_value_type),intent(inout):: this
		if (allocated(this%v)) deallocate(this%v)
	end subroutine dealloc_tarray
	subroutine clean_object(this)
		implicit none
		type(array_value_type):: this
		if (allocated(this%v)) deallocate(this%v)
	end subroutine clean_object
	!
	! construct funktsioonid. extended puhul keerulisem (allokeerimised!)
	subroutine construct_from_value_tarray(this,s,null)
		implicit none
		class(array_value_type),intent(inout):: this
		class(base_value_type),intent(in):: s ! see on allikas
		logical,intent(in),optional:: null ! null value?
		call this%copy(s)
		if (present(null) .and. null) then
			call this%set_null()
		end if
	end subroutine construct_from_value_tarray
	! n - array size
	! ii(1) - esimene indeks, default=1
	subroutine construct_from_params_tarray(this,n,ii,null)
		implicit none
		class(array_value_type),intent(inout):: this
		integer,intent(in):: n
		logical,intent(in),optional:: null ! null value?
		integer,dimension(:),intent(in),optional:: ii ! sellega saab lisamuutujaid sisse anda
		integer:: i1
		! ei kutsu base konstrukti valja, kuna tean, et seal ei tehta midagi!!!
		!call this%base_value_type%construct(n)
		if (allocated(this%v)) deallocate(this%v)
		i1=1
		if (present(ii)) i1=ii(1)
		allocate(this%v(i1:i1+n-1))
		if (present(null) .and. null) then
			call this%set_null()
		end if
	end subroutine construct_from_params_tarray
	!
	! allokeerib sama tyybi... default konstrueerib olemasoleva jargi
	subroutine alloc_same_type_tarray(this,t,construct)
		implicit none
		class(array_value_type),intent(in):: this
		class(base_value_type),intent(inout),allocatable:: t
		logical,intent(in),optional:: construct ! default = .true.
		if (.not.(allocated(t) .and. same_type_as(this,t))) then
			if (allocated(t)) deallocate(t)
			allocate(array_value_type:: t)
		end if
		if (present(construct) .and. .not.construct) return
		call t%construct(this)
	end subroutine alloc_same_type_tarray
	!
	! defineerime vorduse
	! ... saab kasutada ka copy subroutinega
	subroutine copy_value_tarray(this,x)
		implicit none
		class(array_value_type),intent(inout):: this
		class(base_value_type),intent(in):: x
		call this%base_value_type%copy(x)
		select type (x)
		type is (array_value_type)
			if (.not.allocated(x%v)) then
				if (allocated(this%v)) deallocate(this%v)
			else
				if (allocated(this%v)) then
					if (lbound(this%v,dim=1)/=lbound(x%v,dim=1) .or. ubound(this%v,dim=1)/=ubound(x%v,dim=1)) then
						deallocate(this%v)
						allocate(this%v(lbound(x%v,dim=1):ubound(x%v,dim=1)))
					end if
				else
					allocate(this%v(lbound(x%v,dim=1):ubound(x%v,dim=1)))
				end if
				this%v = x%v
			end if
		class default; stop "ERR: copy_value_tarray"
		end select
	end subroutine copy_value_tarray
	subroutine copy_real_tarray(this,x)
		implicit none
		class(array_value_type),intent(inout):: this
		real(rk),intent(in):: x
		call this%base_value_type%copy(x)
		if (allocated(this%v)) then
			this%v = x
			this%null=.false.
		end if
	end subroutine copy_real_tarray
	!
	! defineerime summa vaartuste ja vaartuse ja konstandi vahel
	! ... seda ei ole soovitav kasutada, kuna res allokeeritakse alati uuesti!
	function sum_value_and_value_tarray(this,x) result(res)
		implicit none
		class(array_value_type),intent(in):: this
		class(base_value_type),intent(in):: x
		class(base_value_type),allocatable:: res
		allocate(array_value_type:: res)
		call res%copy(this)
		call res%add(x)
	end function sum_value_and_value_tarray
	function sum_value_and_real_tarray(this,x) result(res)
		implicit none
		class(array_value_type),intent(in):: this
		real(rk),intent(in):: x
		class(base_value_type),allocatable:: res
		allocate(array_value_type:: res)
		call res%copy(this)
		call res%add(x)
	end function sum_value_and_real_tarray
	!
	! kolm voimalust vaartuse lisamiseks olemasolevale vaartusele
	! ... koik kutsutakse valja uhtemoodi
	subroutine add_value_tarray(this,x)
		implicit none
		class(array_value_type),intent(inout):: this
		class(base_value_type),intent(in):: x
		real(rk),dimension(:),allocatable:: y
		integer:: i1,i2
		call this%base_value_type%add(x)
		if (.not.this%null .and. .not.allocated(this%v)) return
		select type (x)
		type is (array_value_type)
			if (.not.allocated(x%v)) return
			if (this%null) then
				if (allocated(this%v)) then
					if (lbound(this%v,dim=1)/=lbound(x%v,dim=1) .or. ubound(this%v,dim=1)/=ubound(x%v,dim=1)) then
						deallocate(this%v)
						allocate(this%v(lbound(x%v,dim=1):ubound(x%v,dim=1)))
					end if
				else
					allocate(this%v(lbound(x%v,dim=1):ubound(x%v,dim=1)))
				end if
				this%v=0.0_rk
			end if
			! vaartus antakse edasi...unset null
			this%null=.false.
			if (lbound(this%v,dim=1)/=lbound(x%v,dim=1) .or. ubound(this%v,dim=1)/=ubound(x%v,dim=1)) then
				i1=min(lbound(this%v,dim=1),lbound(x%v,dim=1))
				i2=max(ubound(this%v,dim=1),ubound(x%v,dim=1))
				allocate(y(lbound(this%v,dim=1):ubound(this%v,dim=1)))
				y=this%v
				deallocate(this%v)
				allocate(this%v(i1:i2))
				this%v=0.0_rk
				this%v(lbound(y,dim=1):ubound(y,dim=1)) = y
				this%v(lbound(x%v,dim=1):ubound(x%v,dim=1)) = this%v(lbound(x%v,dim=1):ubound(x%v,dim=1)) + x%v
				deallocate(y)
			else
				this%v = this%v + x%v
			end if
		class default; stop "ERR: add_value_tarray"
		end select
	end subroutine add_value_tarray
	subroutine add_real_tarray(this,x)
		implicit none
		class(array_value_type),intent(inout):: this
		real(rk),intent(in):: x
		call this%base_value_type%add(x)
		if (allocated(this%v)) then
			this%v=this%v+x
			this%null=.false.
		end if
	end subroutine add_real_tarray
	subroutine add_value_times_constant_tarray(this,x,c)
		implicit none
		class(array_value_type),intent(inout):: this
		class(base_value_type),intent(in):: x
		real(rk),intent(in)::c
		type(array_value_type):: dum
		if (.not.this%null .and. .not.allocated(this%v)) then
			call this%base_value_type%add(x,c)
			return
		end if
		select type (x)
		type is (array_value_type)
			if (.not.allocated(x%v)) then
				call this%base_value_type%add(x,c)
				return
			end if
		class default; stop "ERR: add_value_times_constant_tarray"
		end select
		call dum%copy(x)
		call dum%multiply(c)
		call this%add(dum)
	end subroutine add_value_times_constant_tarray
	!
	! defineerime konstandiga korrutamise
	subroutine multiply_with_constant_tarray(this,c)
		implicit none
		class(array_value_type),intent(inout):: this
		real(rk),intent(in):: c
		call this%base_value_type%multiply(c)
		if (allocated(this%v)) then
			this%v = this%v * c
		end if
	end subroutine multiply_with_constant_tarray
end module array_value_type_module