module base_value_type_module
	use realkind
	implicit none
	private
	public:: base_value_type, accuracy_base_type
	!
	! defineerime funktsiooni vaartuste tapsuse... seda saab hiljem laiendada kui vaja
	type :: accuracy_base_type
		real(rk) :: diff_rel ! see on suhteline viga
		real(rk) :: diff_abs=-1.0 ! see on absoluutne viga (default ei kasutata!)
	end type accuracy_base_type
	!
	type :: base_value_type
		real(rk):: a ! vaartus
		logical:: null=.false. ! vajalik extended puhul... maarab ara kui funktsioon on nullitud!
		! ... null=.true. puhul on naiteks liitmine erinevalt defineeritud osade extended tyypidel
		! ... null vaartust ei kasutada/muudeta (ei tohi kasutada!) base tyybi puhul
	contains
		procedure,pass(this),public:: alloc_same_type => alloc_same_type_tbase
		procedure,pass(this),public:: set_null => set_null_tbase
		procedure,pass(this),public:: dealloc => dealloc_tbase
		procedure,pass(this),public:: is_similar_with => is_similar_with_tbase
		! moned private protseduurid, mis on generic all defineeritud
		procedure,pass(this),private:: construct_from_value => construct_from_value_tbase
		procedure,pass(this),private:: construct_from_params => construct_from_params_tbase
		! protseduurid, mis on vajalikud generic funktsioonide jaoks
		procedure,pass(this),private:: copy_value => copy_value_tbase
		procedure,pass(this),private:: copy_real => copy_real_tbase
		procedure,pass(this),private:: sum_value_and_value => sum_value_and_value_tbase
		procedure,pass(this),private:: sum_value_and_real => sum_value_and_real_tbase
		procedure,pass(this),private:: add_value => add_value_tbase
		procedure,pass(this),private:: add_real => add_real_tbase
		procedure,pass(this),private:: add_value_times_constant => add_value_times_constant_tbase
		procedure,pass(this),private:: multiply_with_constant => multiply_with_constant_tbase
		! define generic procedures
		generic,public:: copy => copy_value, copy_real
		generic,public:: add => add_value, add_real, add_value_times_constant
		generic,public:: multiply => multiply_with_constant
		generic,public:: construct => construct_from_value, construct_from_params
		generic,public:: assignment(=) => copy_value,copy_real
		generic,public:: operator(+) => sum_value_and_value, sum_value_and_real
		! extended tyypide ja allokeerivate muutujate korral
		! ... tuleb defineerida final funktsioon
		!final:: clean
	end type base_value_type
	!
contains
	!
	! define kahe vaartuse tapsuse hindamine
	! isel -- sellega saab erinevaid skeeme implementeerida ning nende vahel valida
	function is_similar_with_tbase(this,ref,acc,isel) result(res)
		implicit none
		real(rk),parameter:: tiny_value=epsilon(1.0)*1.0E3
		class(base_value_type),intent(in):: this
		class(base_value_type),intent(in):: ref
		class(accuracy_base_type),intent(in):: acc
		integer,intent(in),optional:: isel ! sellega saab valida, millist eeskirja kasutatakse!
		logical:: res
		! sanity check... numerical accuracy
		if (abs(this%a)<tiny_value .or. abs(ref%a)<tiny_value .or. abs(this%a - ref%a)<tiny_value) then
			res=.true.
			return
		end if
		!
		! relative difference
		res = abs(this%a - ref%a) <= acc%diff_rel*abs(ref%a)
		!res = abs(this%a - ref%a) <= acc%diff_rel*abs(this%a)
		! absolute difference
		if (acc%diff_abs>0.0_rk) then
			!res = res .or. abs(this%a - ref%a) < acc%diff_abs
			res = res .or. abs(this%a) < acc%diff_abs
		end if
	end function is_similar_with_tbase
	!
	! set null value... paneb ka null parameetri nulliks!
	subroutine set_null_tbase(this)
		implicit none
		class(base_value_type),intent(inout):: this
		this%null=.true.
		this%a=0.0_rk
	end subroutine set_null_tbase
	! see on ainult kohahoidja
	subroutine dealloc_tbase(this)
		implicit none
		class(base_value_type),intent(inout):: this
	end subroutine dealloc_tbase
	!
	! construct funktsioonid. extended puhul keerulisem (allokeerimised!)
	subroutine construct_from_value_tbase(this,s,null)
		implicit none
		class(base_value_type),intent(inout):: this
		class(base_value_type),intent(in):: s ! see on allikas
		logical,intent(in),optional:: null ! null value?
		if (present(null) .and. null) then
			call this%set_null()
		else
			call this%copy(s)
		end if
	end subroutine construct_from_value_tbase
	subroutine construct_from_params_tbase(this,n,ii,null)
		implicit none
		class(base_value_type),intent(inout):: this
		integer,intent(in):: n
		logical,intent(in),optional:: null ! null value?
		integer,dimension(:),intent(in),optional:: ii ! sellega saab lisamuutujaid sisse anda
		if (present(null) .and. null) then
			call this%set_null()
		end if
	end subroutine construct_from_params_tbase
	!
	! allokeerib sama tyybi... default konstrueerib olemasoleva jargi
	subroutine alloc_same_type_tbase(this,t,construct)
		implicit none
		class(base_value_type),intent(in):: this
		class(base_value_type),intent(inout),allocatable:: t
		logical,intent(in),optional:: construct
		if (.not.(allocated(t) .and. same_type_as(this,t))) then
			if (allocated(t)) deallocate(t)
			allocate(base_value_type:: t)
		end if
		if (present(construct) .and. .not.construct) return
		call t%construct(this)
	end subroutine alloc_same_type_tbase
	!
	! defineerime vorduse
	! ... saab kasutada ka copy subroutinega
	subroutine copy_value_tbase(this,x)
		implicit none
		class(base_value_type),intent(inout):: this
		class(base_value_type),intent(in):: x
		this%a=x%a
		this%null=x%null
	end subroutine copy_value_tbase
	subroutine copy_real_tbase(this,x)
		implicit none
		class(base_value_type),intent(inout):: this
		real(rk),intent(in):: x
		this%a=x
	end subroutine copy_real_tbase
	!
	! defineerime summa vaartuste ja vaartuse ja konstandi vahel
	! ... seda ei ole soovitav kasutada, kuna res allokeeritakse alati uuesti!
	function sum_value_and_value_tbase(this,x) result(res)
		implicit none
		class(base_value_type),intent(in):: this
		class(base_value_type),intent(in):: x
		class(base_value_type),allocatable:: res
		res%a=this%a + x%a
	end function sum_value_and_value_tbase
	function sum_value_and_real_tbase(this,x) result(res)
		implicit none
		class(base_value_type),intent(in):: this
		real(rk),intent(in):: x
		class(base_value_type),allocatable:: res
		res%a=this%a + x
	end function sum_value_and_real_tbase
	!
	! kolm voimalust vaartuse lisamiseks olemasolevale vaartusele
	! ... koik kutsutakse valja uhtemoodi
	subroutine add_value_tbase(this,x)
		implicit none
		class(base_value_type),intent(inout):: this
		class(base_value_type),intent(in):: x
		this%a=this%a+x%a
	end subroutine add_value_tbase
	subroutine add_real_tbase(this,x)
		implicit none
		class(base_value_type),intent(inout):: this
		real(rk),intent(in):: x
		this%a=this%a+x
	end subroutine add_real_tbase
	subroutine add_value_times_constant_tbase(this,x,c)
		implicit none
		class(base_value_type),intent(inout):: this
		class(base_value_type),intent(in):: x
		real(rk),intent(in)::c
		this%a=this%a+x%a*c
	end subroutine add_value_times_constant_tbase
	!
	! defineerime konstandiga korrutamise
	subroutine multiply_with_constant_tbase(this,c)
		implicit none
		class(base_value_type),intent(inout):: this
		real(rk),intent(in):: c
		this%a=this%a*c
	end subroutine multiply_with_constant_tbase
end module base_value_type_module