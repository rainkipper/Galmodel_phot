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
		
		IS = .false.
		mmodal = .false. 
		nlive = 2 !testiks nii v2ike
		ceff = .true.
		tol = 0.5 !ei tea, mis siia peaks k2ima
		efr = 0.8
		ndims = leia_Ndim()!hiljem
		nPar = ndims !hiljem
		nCdims = 1
		maxModes = 1
		updInt  = 1
		Ztol = -1.d90
		root = "Output/"
		seed = 12345
		allocate(pWrap(1:nPar)); pWrap = 0 !
		feedback = .false.
		resume = .false.
		outfile = .true.
		initMPI = .false.
		logZero = -1.0d10 !ei ole kindel selles
		maxiter = 1 !ehk umbes minut
		context = 0 !mittevajalik
		
		
		!
		! ========== asjade initsialiseerimine... sh all_comp jm ===============
		!
		call convert_input_comp_to_all_comp(input_comps, all_comp)
		call asenda_viited(input_comps, all_comp)
		
		!
		! =============== fittimine ise ================
		!

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
			double precision Cube(n_dim), lnew
			integer :: i
			type(prof_par_list_type), pointer ::  par_list
			integer mitmes_cube
			! all_comp uuendamine
			! algselt uuendab input_comp v22rtused ning siis edasi all_comp
			mitmes_cube = 0
			do i=1,size(input_comps)
				if(input_comps(i)%incl%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%incl%val = UniformPrior(cube(mitmes_cube), input_comps(i)%incl%min,  input_comps(i)%incl%max)
				end if	
				if(input_comps(i)%cnt_x%kas_fitib) then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%cnt_x%val = UniformPrior(cube(mitmes_cube), input_comps(i)%cnt_x%min,  input_comps(i)%cnt_x%max)
				end if	
				if(input_comps(i)%cnt_y%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%cnt_y%val = UniformPrior(cube(mitmes_cube), input_comps(i)%cnt_y%min,  input_comps(i)%cnt_y%max)
				end if	
				if(input_comps(i)%pos%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%pos%val = UniformPrior(cube(mitmes_cube), input_comps(i)%pos%min,  input_comps(i)%pos%max)
				end if	
				if(input_comps(i)%theta0%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					input_comps(i)%theta0%val = UniformPrior(cube(mitmes_cube), input_comps(i)%theta0%min,  input_comps(i)%theta0%max)
				end if	
				par_list=>input_comps(i)%prof_pars
				do while(par_list%filled)
					if(par_list%par%kas_fitib)  then
					mitmes_cube=mitmes_cube+1
					par_list%par%val = UniformPrior(cube(mitmes_cube), par_list%par%min,  par_list%par%max)
				end if	
					if(associated(par_list%next)) then
						par_list => par_list%next
					else
						exit
					end if
				end do
			end do

			call asenda_viited(input_comps, all_comp) !all_comp muutujas asendamine


			!loglike reaalne arvutamine
			lnew = calc_log_likelihood(all_comp, images)

		end subroutine fun_loglike
		subroutine fun_dumper(nSamples,nlive,nPar,physLive,posterior, paramConstr,maxloglike,logZ,INSlogZ,logZerr,context)
			implicit none
			integer :: nlive, nSamples, nPar, context
			double precision :: maxloglike, logZ, INSlogZ,logZerr
			double precision, pointer :: posterior(:,:)
			double precision, pointer :: physLive(:,:)
			double precision, pointer :: paramConstr(:)
! 			double precision, dimension(1:nlive, 1:nPar+1) :: physLive
! 			double precision, dimension(1, 1:4*nPar) :: paramConstr
			print*, "Ikka elus..."
		end subroutine fun_dumper
end subroutine jooksuta_fittimine
	
	!		prior fyysikalisesse vahemikku
	! 		function UniformPrior(r,x1,x2)
	! 		      	implicit none
	! 		      	double precision r,x1,x2,UniformPrior
	! 		      	UniformPrior=x1+r*(x2-x1)
	! 		end function UniformPrior
end module fitting_multinest_module
