module fitting_multinest_module
	use nested
	use likelihood_module
	use output_module
	real(rk), dimension(:), allocatable :: parim_phys_points !ehk outputi jaoks j2tab viimase siia
contains
	subroutine fittimine_multinest(images, input_comps, all_comp)
		implicit none
		!
		! jagatud sisend-v2ljund teiste fittijatega
		!
		type(comp_input_type), dimension(:), allocatable, intent(inout), target :: input_comps
		type(all_comp_type), intent(inout) :: all_comp !ehk v2ljundiks
		type(image_type), dimension(:), allocatable, intent(in) :: images
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
! 		logical, dimension(:), allocatable :: recalc_comp
		
		IS = .true.
		mmodal = .false. 
		ceff = .true.
		tol = 0.5
		efr = multinest_efr
		ndims = leia_vabade_parameetrite_arv(input_comps) !moodulist comp.f90
		nlive = ndims + N_multinest_extra_points !testiks nii v2ike
		nPar = ndims !hiljem
		nCdims = 1
		maxModes = 3 !rohkem on mottetu, kui tahta hiljem multinesti rakendada
		updInt  = 3
		Ztol = -1.d90
		root = repeat(" ",100); root = trim(multinest_output_header) !m2lu korrastamiseks alguses
		seed = -1
		allocate(pWrap(1:nPar)); pWrap = 0 !
		feedback = .false.
		resume = .false.
		outfile = .true.
		initMPI = .false.
		logZero = -1.0d10 !ei ole kindel selles
		maxiter = 2000 !et saada ainult esialgne l2hend
		context = 0 !mittevajalik
		
		
		!
		! ========== asjade initsialiseerimine... sh all_comp jm ===============
		!
! 		allocate(recalc_comp(1:all_comp%N_comp)); recalc_comp = .true. !ehk koik peab esimene kord uuesti arvutama... likelihoodi jaoks on see globaalne parameeter, mida pidevalt muudetakse

		
		!
		! =============== fittimine ise ================
		!
		
		
		call nestRun(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nCdims, maxModes, updInt, Ztol, root, seed, &
			 pWrap, feedback, resume, outfile, initMPI, logZero, maxiter, fun_loglike, fun_dumper, context)

		
	contains

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
			nullify(par_list)

			call convert_input_comp_to_all_comp(input_comps, all_comp)
			call asenda_viited(input_comps, all_comp) !all_comp muutujas asendamine
			lnew =  calc_log_likelihood(all_comp, images)
			nullify(par_list)
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
			if(.not.allocated(parim_phys_points)) allocate(parim_phys_points(1:nPar))
			parim_phys_points = physLive(parim, 1:nPar) !salvestab parimad punktid mooduli muutujaks, et hiljem output oleks lihtsam
			if(.not.kas_vaikselt) print*, "============== parim on ", maxloglike, size(physLive)
			if(.not.kas_vaikselt) print*, "(A,2F15.6)", "============== mean, sd of LL", keskmine, sd_h2lve
			mitmes_cube = 0
			do i=1,size(input_comps)
				if(input_comps(i)%incl%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%incl%val = physLive(parim,mitmes_cube)
				end if
				if(input_comps(i)%cnt_x%kas_fitib) then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%cnt_x%val =  physLive(parim,mitmes_cube)
				end if
				if(input_comps(i)%cnt_y%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%cnt_y%val =  physLive(parim,mitmes_cube)
				end if
				if(input_comps(i)%pos%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%pos%val = physLive(parim,mitmes_cube)
				end if
				if(input_comps(i)%theta0%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%theta0%val = physLive(parim,mitmes_cube)
				end if
				par_list=>input_comps(i)%prof_pars
				do while(par_list%filled)
					if(par_list%par%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					par_list%par%val = physLive(parim,mitmes_cube)
				end if
					if(associated(par_list%next)) then
						par_list => par_list%next
					else
						exit
					end if
				end do
			end do
			nullify(par_list)
			call output_images(input_comps, images, all_comp)
			if(mis_fittimise_tyyp == 2) call output_ML(input_comps, images)
			call output_like_input(input_comps)
			if(.not.kas_vaikselt) print "(A,F15.10)", "========================================================== t per LL = ", (dt-alguse_aeg)/LL_counter
		end subroutine fun_dumper
	end subroutine fittimine_multinest
		
end module fitting_multinest_module
