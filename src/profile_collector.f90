include "src/Profiles/prof_Einasto.f90"
include "src/Profiles/prof_Sersic.f90"
include "src/Profiles/prof_dustplane.f90" 
module profile_collector_module
	use prof_Einasto_module
	use prof_Sersic_module
	use prof_dustplane_module
end module profile_collector_module