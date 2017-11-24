module comp_image_module	
	use pixel_module
	use images_module
	use yldine_matemaatika_module
	!peab kasutama veel pildi infot
	!AINULT ruudulisele gridile
	
	!nurkade koordinaadid
	real(rk), parameter, dimension(1:4), private :: pix_nihe_x = [-1.0, -1.0, 0.0, 0.0] !m22rab pikslite j2rjestuse massiivis (vt pixel_module)
	real(rk), parameter, dimension(1:4), private :: pix_nihe_y = [-1.0, 0.0, 0.0, -1.0] 
! 	real(rk), parameter, dimension(1:4), private :: pix_nihe_x = [-0.5, -0.5, 0.5, 0.5] !kui piksli v22rtus keskel voetud
! 	real(rk), parameter, dimension(1:4), private :: pix_nihe_y = [-0.5, 0.5, 0.5, -0.5] !kui piksli v22rtus keskel voetud
! 	real(rk), parameter, dimension(1:4), private :: pix_nihe_x = [0.0, 0.0, 1.0, 1.0]
! 	real(rk), parameter, dimension(1:4), private :: pix_nihe_y = [0.0, 1.0, 1.0, 0.0]

	
	type :: comp_image_real_type 
		type(square_pixel_type), dimension(:,:), allocatable :: pix
		real(rk), dimension(:,:), allocatable :: mx !reaalne maatriks tuleb siia
		real(rk), dimension(:,:), allocatable :: M_enne_tasandit, M_p2rast_tasandit !enne ja p2rast tolmu tasandit piki vaatejoont
		!
		logical :: recalc_XcYc_coords !kas vaja arvutada fyysikalistest koordinaatidest uuesti komponenti koordinaadid
		logical :: recalc_image
		integer :: corresponding_adaptive_im, corresponding_adaptive_im_enne_tasandit, corresponding_adaptive_im_p2rast_tasandit	!kui kasutab adaptive_im abi, et pilti joonistada
	end type comp_image_real_type
contains
	subroutine create_comp_image_from_obs(mdl, im)
		!teeb raamistiku vaatluspildi pohjal ning lisab koordinaadid ja algsed teisendused
		implicit none
		type(comp_image_real_type), intent(out) :: mdl
		type(image_type), intent(in) :: im
		integer :: Nx, Ny, i,j

		Nx = size(im%obs, 1); Ny = size(im%obs, 2)
		if(.not.allocated(mdl%pix)) allocate(mdl%pix(1:Nx, 1:Ny))
		if(.not.allocated(mdl%mx)) allocate(mdl%mx(1:Nx, 1:Ny))
		if(.not.allocated(mdl%M_enne_tasandit)) allocate(mdl%M_enne_tasandit(1:Nx, 1:Ny))
		if(.not.allocated(mdl%M_p2rast_tasandit)) allocate(mdl%M_p2rast_tasandit(1:Nx, 1:Ny))
		mdl%recalc_image = .true. !ehk see koige algsem ning midagi pole veel arvutatud
		mdl%recalc_XcYc_coords = .true.
		mdl%corresponding_adaptive_im = -1 !ehk pole

		!vaja teha ainult sisselugemise aeg, seega kiirus pole oluline
		do i=1,Nx
		do j=1,Ny
			mdl%pix(i,j)%Xi_nurgad = pix_nihe_x+real(i, kind=rk)
			mdl%pix(i,j)%Yi_nurgad = pix_nihe_y+real(j, kind=rk)
			call im%XiYi_to_XpYp(mdl%pix(i,j)%Xi_nurgad, mdl%pix(i,j)%Yi_nurgad, mdl%pix(i,j)%Xp_nurgad, mdl%pix(i,j)%Yp_nurgad)
		end do
		end do
	end subroutine create_comp_image_from_obs
	
	

	subroutine fill_corners(mdl, func, func_approx) !t2itab func v22rtusega
		!vaja, et pikslite komponendi koordinaadid oleks juba oiged
		implicit none
		interface
			function func(Xc, Yc) result(res)
				import rk
				real(rk), intent(in) :: Xc, Yc
				real(rk) :: res
			end function func
		end interface
		interface
			function func_approx(Xc, Yc) result(res)
				import rk
				real(rk), intent(in) :: Xc, Yc
				real(rk) :: res
			end function func_approx
		end interface
		type(comp_image_real_type), intent(inout) :: mdl
		integer :: i,j, nurk
		integer :: Nx, Ny
		real(rk) :: v22rtus
		Nx = size(mdl%pix, 1)
		Ny = size(mdl%pix, 2)
		!
		! ========== koordinaatide t2itmine =============
		!
		!keskosa t2itmine
		do i = 2,Nx
		do j = 2,Ny
			if(mdl%pix(i,j)%Xc_nurgad(1)**2 + mdl%pix(i,j)%Yc_nurgad(1) **2 < adaptive_image_dist_piirang**2) then
				v22rtus = func( mdl%pix(i,j)%Xc_nurgad(1), mdl%pix(i,j)%Yc_nurgad(1) )
			else
				v22rtus = func_approx( mdl%pix(i,j)%Xc_nurgad(1), mdl%pix(i,j)%Yc_nurgad(1) )
			end if
			mdl%pix(i, j)%val_nurgad(1) = v22rtus
			mdl%pix(i-1, j-1)%val_nurgad(3) = v22rtus
			mdl%pix(i-1, j)%val_nurgad(4) = v22rtus
			mdl%pix(i, j-1)%val_nurgad(2) = v22rtus
		end do
		end do
		!servad k2sitsi yle teha
		do i=1,Nx
			do nurk = 1, 4
				mdl%pix(i, 1)%val_nurgad(nurk) = func_approx(mdl%pix(i, 1)%Xc_nurgad(nurk), mdl%pix(i, 1)%Yc_nurgad(nurk))
				mdl%pix(i, Ny)%val_nurgad(nurk) = func_approx(mdl%pix(i, Ny)%Xc_nurgad(nurk), mdl%pix(i, Ny)%Yc_nurgad(nurk))
			end do
		end do
		do j=1,Ny
			do nurk = 1,4
				mdl%pix(1, j)%val_nurgad(nurk) = func_approx(mdl%pix(1, j)%Xc_nurgad(nurk), mdl%pix(1, j)%Yc_nurgad(nurk))
				mdl%pix(Nx, j)%val_nurgad(nurk) = func_approx(mdl%pix(Nx, j)%Xc_nurgad(nurk), mdl%pix(Nx, j)%Yc_nurgad(nurk))
			end do
		end do
	end subroutine fill_corners
	subroutine fill_corners_exact(mdl, func, func_approx)
		implicit none
		interface
			function func(Xc, Yc) result(res)
				import rk
				real(rk), intent(in) :: Xc, Yc
				real(rk) :: res
			end function func
			function func_approx(Xc, Yc) result(res)
				import rk
				real(rk), intent(in) :: Xc, Yc
				real(rk) :: res
			end function func_approx
		end interface
		type(comp_image_real_type), intent(inout) :: mdl
		integer :: i,j
		
		do i = 1,size(mdl%pix,1)
		do j = 1,size(mdl%pix,2)
			if(kas_kasutab_l2hendit_piksli_v22rtuse_jaoks(mdl%pix(i,j))) then
				mdl%pix(i,j)%val_nurgad(1) = func_approx(mdl%pix(i,j)%Xc_nurgad(1), mdl%pix(i,j)%Yc_nurgad(1))
				mdl%pix(i,j)%val_nurgad(2) = func_approx(mdl%pix(i,j)%Xc_nurgad(2), mdl%pix(i,j)%Yc_nurgad(2))
				mdl%pix(i,j)%val_nurgad(3) = func_approx(mdl%pix(i,j)%Xc_nurgad(3), mdl%pix(i,j)%Yc_nurgad(3))
				mdl%pix(i,j)%val_nurgad(4) = func_approx(mdl%pix(i,j)%Xc_nurgad(4), mdl%pix(i,j)%Yc_nurgad(4))
			else
				mdl%pix(i,j)%val_nurgad(1) = func(mdl%pix(i,j)%Xc_nurgad(1), mdl%pix(i,j)%Yc_nurgad(1))
				mdl%pix(i,j)%val_nurgad(2) = func(mdl%pix(i,j)%Xc_nurgad(2), mdl%pix(i,j)%Yc_nurgad(2))
				mdl%pix(i,j)%val_nurgad(3) = func(mdl%pix(i,j)%Xc_nurgad(3), mdl%pix(i,j)%Yc_nurgad(3))
				mdl%pix(i,j)%val_nurgad(4) = func(mdl%pix(i,j)%Xc_nurgad(4), mdl%pix(i,j)%Yc_nurgad(4))				
			end if
		end do
		end do
		
	end subroutine fill_corners_exact
	
	function combine_comp_images_to_make_image(mdl, weights) result(res)
		implicit none
		type(comp_image_real_type), intent(in), dimension(:), allocatable :: mdl
		real(rk), dimension(:,:), allocatable :: res
		real(rk), dimension(:), allocatable, intent(in) :: weights !praktikas on need mass-heledus suhted
		integer :: i
		!andmete kontroll
		if(size(weights,1) .ne. size(mdl, 1)) then
			print*, "Kaalude ja mudelpiltide massiivi pikkus peab sama olema"
			stop
		end if
		
		!funktsioon ise on siit...
		allocate(res(1:size(mdl(1)%pix, 1), 1:size(mdl(1)%pix, 2)))
		res = 0.0_rk
		do i=1,size(mdl, 1)
			!kontroll, kas pikslite suurused on ikkagi samad igal juhul
			if(size(mdl(i)%pix, 1)==size(res,1) .and. size(mdl(i)%pix, 2)==size(res,2)) then 
				res = res + mdl(i)%mx * weights(i)
			else
				print*, "Pildi suurused ei klapi mudelpiltide kokkupanekul"
				stop
			end if
		end do
	end function combine_comp_images_to_make_image
	
end module comp_image_module
	
	