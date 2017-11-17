module populations_module
	use yldine_matemaatika_module
	use file_operations_module
	type :: population_type
		type(distribution_1D_type) :: spec
		character(default_character_length) :: name
	end type population_type
	type(population_type), dimension(:), allocatable, target :: populations

	
contains
	subroutine init_populations(populations_file)
		implicit none
		character(len=default_character_length), intent(in) :: populations_file
		
		call init_populations_Blanton(populations_file)
		
	end subroutine init_populations
	
	subroutine init_populations_Blanton(populations_file)
		implicit none
		real(rk), dimension(:,:), allocatable :: mx
		character(len=default_character_length), intent(in) :: populations_file
		character(len=default_character_length), parameter, dimension(1:5) :: nimed = ["Blanton_1", "Blanton_2", "Blanton_3", "Blanton_4", "Blanton_5"]
		integer :: i
		mx = read_tabel(populations_file, size(nimed,1)+1)
		allocate(populations(1:size(nimed,1)))
		do i=1,size(nimed, 1)
			populations(i)%name = nimed(i)
			allocate(populations(i)%spec%x(1:size(mx,1)))
			allocate(populations(i)%spec%f(1:size(mx,1)))
			populations(i)%spec%x = mx(:,1)
			populations(i)%spec%f = mx(:,i+1)
		end do
	end subroutine init_populations_Blanton
end  module populations_module