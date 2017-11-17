module dust_module
	use yldine_matemaatika_module
	!tolmu parameetrid voetud Calzetti juurest http://iopscience.iop.org/article/10.1086/324269/pdf valemid 8a ja 8b
contains
	elemental function dust_kappa(l) result(res)
		implicit none
		real(rk), intent(in) :: l
		real(rk) :: mu, res
		mu = l * 1e-4 !valemid on antud mikronites... vaja angstromides
		if(mu<0.12 .or. mu>2.20) then
			res = 0
		else
			if(mu < 0.63) then
				res = 1.17*(-2.156 + 1.509/mu - 0.198/(mu*mu) + 0.011/(mu**3)) + 1.78
			else
				res = 1.17*(-1.857 + 1.040/mu) + 1.78
			end if
		end if
	end function dust_kappa
end module dust_module