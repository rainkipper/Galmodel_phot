module fitting_multinest_module
	use nested
	use likelihood_module

contains
	subroutine jooksuta_fittimine(images, input_comps, all_comp)
		implicit none
		!
		! jagatud sisend-v2ljund teiste fittijatega
		!
		type(comp_input_type), dimension(:), allocatable, intent(inout), target :: input_comps
		type(all_comp_type), intent(out) :: all_comp !ehk v2ljundiks
		type(image_type), dimension(:), allocatable, intent(in) :: images
			real(rk) :: alguse_aeg
		!
		! ======= multinesti muutujad =======
		!
		logical :: IS
		logical :: mmodal
		integer :: nlive
		logical :: ceff
		double precision :: tol
		double precision :: efr
		integer :: ndims
		integer :: nPar
		integer  :: nCdims
		integer  :: maxModes
		integer  :: updInt
		double precision :: Ztol
		character(LEN=100)  :: root
		integer  :: seed
		integer, dimension(:), allocatable  :: pWrap
		logical  :: feedback
		logical  :: resume
		logical  :: outfile
		logical  :: initMPI
! 		double precision  :: 
		integer  :: maxiter
		integer ::  context
		double precision :: logZero
		logical, dimension(:), allocatable :: recalc_comp
		
		IS = .true.
		mmodal = .false. 
		nlive = 15 !testiks nii v2ike
		ceff = .true.
		tol = 0.5
		efr = 0.99
		ndims = leia_Ndim()!hiljem
		nPar = ndims !hiljem
		nCdims = 1
		maxModes = 1
		updInt  = 1
		Ztol = -1.d90
		root = "Output/"
		seed = -1
		allocate(pWrap(1:nPar)); pWrap = 0 !
		feedback = .false.
		resume = .false.
		outfile = .true.
		initMPI = .false.
		logZero = -1.0d10 !ei ole kindel selles
		maxiter = 10000
		context = 0 !mittevajalik
		
		
		!
		! ========== asjade initsialiseerimine... sh all_comp jm ===============
		!
		allocate(recalc_comp(1:all_comp%N_comp)); recalc_comp = .true. !ehk koik peab esimene kord uuesti arvutama... likelihoodi jaoks on see globaalne parameeter, mida pidevalt muudetakse
		call convert_input_comp_to_all_comp(input_comps, all_comp)
		call asenda_viited(input_comps, all_comp)
		all_comp%comp(:)%adaptive_image_number = -1 !default -1, et teeks uue adaptiivse pildi esimene kord
		call init_calc_log_likelihood(all_comp, images) !s2ttib likelihoodi mooduli muutujad, et v2hendada arvutamisis
		
		!
		! =============== fittimine ise ================
		!
		
		call cpu_time(alguse_aeg)
		call nestRun(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nCdims, maxModes, updInt, Ztol, root, seed, &
			 pWrap, feedback, resume, outfile, initMPI, logZero, maxiter, fun_loglike, fun_dumper, context)
		
	contains
		function leia_Ndim() result(res)
			implicit none
			integer :: i
			type(prof_par_list_type), pointer ::  par_list
			integer :: res
			res = 0
			do i=1,size(input_comps)
				if(input_comps(i)%incl%kas_fitib) res=res+1
				if(input_comps(i)%cnt_x%kas_fitib) res=res+1
				if(input_comps(i)%cnt_y%kas_fitib) res=res+1
				if(input_comps(i)%pos%kas_fitib) res=res+1
				if(input_comps(i)%theta0%kas_fitib) res=res+1
				par_list=>input_comps(i)%prof_pars
				do while(par_list%filled)
					if(par_list%par%kas_fitib) res=res+1
					if(associated(par_list%next)) then
						par_list => par_list%next
					else
						exit
					end if
				end do
			end do
			print*, "Ndim = ", res
		end function leia_Ndim
		subroutine fun_loglike(Cube,n_dim,nPar,lnew,context)
			implicit none
			integer ::  n_dim, nPar, context
			double precision :: Cube(n_dim),lnew
			real(rk) :: tmp(4) !ajutisteks test
			integer :: i
			type(prof_par_list_type), pointer ::  par_list
			integer mitmes_cube
			real(rk), dimension(:), allocatable :: lisakaalud_massile
			
			mitmes_cube = 0
			do i=1,size(input_comps)
				if(input_comps(i)%incl%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					cube(mitmes_cube) = UniformPrior(cube(mitmes_cube), input_comps(i)%incl%min,  input_comps(i)%incl%max)
					input_comps(i)%incl%val = cube(mitmes_cube)
				end if	
				if(input_comps(i)%cnt_x%kas_fitib) then
					mitmes_cube=mitmes_cube+1
					cube(mitmes_cube) = UniformPrior(cube(mitmes_cube), input_comps(i)%cnt_x%min,  input_comps(i)%cnt_x%max)
					input_comps(i)%cnt_x%val = cube(mitmes_cube)	
				end if	
				if(input_comps(i)%cnt_y%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					cube(mitmes_cube) = UniformPrior(cube(mitmes_cube), input_comps(i)%cnt_y%min,  input_comps(i)%cnt_y%max)
					input_comps(i)%cnt_y%val = cube(mitmes_cube)
				end if	
				if(input_comps(i)%pos%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					cube(mitmes_cube) = UniformPrior(cube(mitmes_cube), input_comps(i)%pos%min,  input_comps(i)%pos%max)
					input_comps(i)%pos%val = cube(mitmes_cube)
				end if	
				if(input_comps(i)%theta0%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					cube(mitmes_cube) = UniformPrior(cube(mitmes_cube), input_comps(i)%theta0%min,  input_comps(i)%theta0%max)
					input_comps(i)%theta0%val = cube(mitmes_cube)
				end if	
				par_list=>input_comps(i)%prof_pars
				do while(par_list%filled)
					if(par_list%par%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					cube(mitmes_cube) = UniformPrior(cube(mitmes_cube), par_list%par%min,  par_list%par%max)
					par_list%par%val = cube(mitmes_cube)
					
				end if	
					if(associated(par_list%next)) then
						par_list => par_list%next
					else
						exit
					end if
				end do
			end do
			

			call convert_input_comp_to_all_comp(input_comps, all_comp)
			call asenda_viited(input_comps, all_comp) !all_comp muutujas asendamine
			lnew =  calc_log_likelihood(all_comp, images, lisakaalud_massile)
			
			!kui lisamassidele juurde asju arvutatud, siis paneb uued massid vastavalt eelmistele... input_comps juurde
			do i=1,all_comp%N_comp
				par_list => input_comps(i)%prof_pars
				do while(par_list%filled)
					if(trim(par_list%par_name)=="M") then
						par_list%par%val = par_list%par%val * lisakaalud_massile(i)
						exit
					end if
					par_list => par_list%next
				end do
			end do
		end subroutine fun_loglike
		subroutine fun_dumper(nSamples,nlive,nPar,physLive, posterior, paramConstr,maxloglike,logZ,INSlogZ,logZerr,context)
			implicit none
			integer :: nlive, nSamples, nPar, context
			double precision :: maxloglike, logZ, INSlogZ,logZerr
			double precision, pointer :: posterior(:,:)
			double precision, pointer :: physLive(:,:)
			double precision, pointer :: paramConstr(:)
			double precision :: keskmine, sd_h2lve
			integer :: i
			type(prof_par_list_type), pointer ::  par_list
			integer mitmes_cube	
			integer :: parim	
			real(rk) :: dt	
			parim = maxloc(physLive(:,nPar+1),1)
			keskmine = sum(posterior(:,nPar+1))/size(posterior, 1)
			sd_h2lve = sqrt(sum( (posterior(:,nPar+1)-keskmine)**2 )/size(posterior, 1))
			print*, size(posterior, 1), size(posterior, 2)
			print*, "============== parim on ", maxloglike, size(physLive)
			print "(A,2F15.6)", "============== mean, sd of LL", keskmine, sd_h2lve
			mitmes_cube = 0
			do i=1,size(input_comps)
				if(input_comps(i)%incl%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					print*, trim(all_comp%comp(i)%comp_name)," Incl", physLive(parim,mitmes_cube)*180.0/pi
				end if
				if(input_comps(i)%cnt_x%kas_fitib) then
					mitmes_cube=mitmes_cube+1
					print*, trim(all_comp%comp(i)%comp_name)," cnt_x", physLive(parim,mitmes_cube)/arcsec_to_rad
				end if
				if(input_comps(i)%cnt_y%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					print*, trim(all_comp%comp(i)%comp_name)," cnt_y", physLive(parim,mitmes_cube)
				end if
				if(input_comps(i)%pos%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					print*, trim(all_comp%comp(i)%comp_name)," pos", physLive(parim,mitmes_cube)*180/pi
				end if
				if(input_comps(i)%theta0%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					print*, trim(all_comp%comp(i)%comp_name)," theta0",	physLive(parim,mitmes_cube)

				end if
				par_list=>input_comps(i)%prof_pars
				do while(par_list%filled)
					if(par_list%par%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					print*, trim(all_comp%comp(i)%comp_name)," ",trim(par_list%par_name), physLive(parim,mitmes_cube)
				end if
					if(associated(par_list%next)) then
						par_list => par_list%next
					else
						exit
					end if
				end do
			end do
			print*, "-------"
			print*, "Massid:"
			do i=1,size(input_comps)
				par_list=>input_comps(i)%prof_pars
				do while(par_list%filled)
					if(trim(par_list%par_name)=="M") then
						print*, " ", trim(input_comps(i)%comp_name), par_list%par%val
						exit
					end if
					par_list=>par_list%next
				end do
			end do
			call cpu_time(dt)
			print "(A,F15.10)", "========================================================== dt_algusest = ", dt-alguse_aeg
		end subroutine fun_dumper
end subroutine jooksuta_fittimine
	
	!		prior fyysikalisesse vahemikku
	! 		function UniformPrior(r,x1,x2)
	! 		      	implicit none
	! 		      	double precision r,x1,x2,UniformPrior
	! 		      	UniformPrior=x1+r*(x2-x1)
	! 		end function UniformPrior
end module fitting_multinest_module
