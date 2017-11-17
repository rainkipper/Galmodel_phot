module distributions_module
	use constants_module
	type :: distribution_1D_type
		sequence
		real(rk), dimension(:), allocatable :: x 
		real(rk), dimension(:), allocatable :: f
	end type distribution_1D_type
contains


end module distributions_module