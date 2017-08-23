module comp_image_real_module
	use yldine_matemaatika_module !vaja interpoleerimiseks
	use constants_module
	implicit none
	
	integer, parameter :: im_list_maxsize = 101
	integer, parameter :: puudub = -199
	integer, parameter :: maxsize = 30000 !lambi suurus, et mahutada punktide massiive
	real(rk), parameter :: x0_default = -35.0_rk
	real(rk), parameter :: y0_default = -35.0_rk
	real(rk), parameter :: x1_default =  35.0_rk
	real(rk), parameter :: y1_default =  35.0_rk
	real(rk), private :: comp_im_edasijagamise_maksimaalne_abs_t2psus = 1.000001_rk !peab hiljem automaatselt t2itma
	integer, private :: comp_im_maxlevel = 12
	real(rk), parameter :: comp_im_edasijagamise_threshold = 0.04
! 	integer :: comp_im_kokku = 1
	

	type :: rk_point_type
		real(rk), pointer :: point
	end type rk_point_type
	
	type :: logical_point_type !vaja lao massiivide jaoks... fortrani erip2ra
		logical, pointer :: point
	end type logical_point_type
	
	type :: comp_image_type
		real(rk), dimension(1:5) 				:: x,y !jrk on bottom-left, top-left, top-right, bottom-right, centre
		type(rk_point_type), dimension(1:5)	:: val !v22rtus funktsioonile... voib olla suvaline asi
		integer, dimension(1:5)					:: id
! 		real(rk)					:: sum_val !summaarne heledus kasti sees
! 		real(rk)					:: sum_val_korda_pind !reaalne v22rtus, mis voetakse kui koik subpixslis
		logical 								:: last_level = .true. !kas omadega pohjas
		logical 								:: kas_paremvasak !kas subkomponent on paremal-vasak poolitatud voi yles/alla
		type(comp_image_type), pointer 			:: sub1, sub2 !viited alamstruktuurile
		contains
			procedure :: get_val => get_comp_im_val
	end type comp_image_type
	
	type :: comp_image_linked_list_type
		integer :: mitmes_see_comp_im_on = 1
		type(comp_image_type)							:: comp_im
		real(rk), dimension(1:maxsize)		:: val_ladu
		real(rk), dimension(1:maxsize)					:: x_ladu, y_ladu
		logical, dimension(1:maxsize) 					:: kas_arvutatud_ladu
		integer, dimension(1:maxsize)					:: id_list
		type(comp_image_linked_list_type), pointer 		:: next
	end type comp_image_linked_list_type
	
	type(comp_image_linked_list_type), target 	:: comp_im_list !globaalne muutuja, kuhu koik sisse tulevad
	integer 							:: comp_im_total_counter !peab arvet, palju on kokku muutujaid... saaks ka lugedes alati
	

	

	
	contains
		
		recursive function get_comp_im_val(comp_im, Xc, Yc) result(res)
			implicit none
			class(comp_image_type), intent(in) :: comp_im
			real(rk), intent(in) :: Xc, Yc
			real(rk) :: res
			real(rk) :: f_low, f_up
			if(comp_im%last_level) then
				f_low = (comp_im%val(4)%point - comp_im%val(1)%point)/(comp_im%x(4) - comp_im%x(1)) * (Xc-comp_im%x(1)) + comp_im%val(1)%point
				f_up = (comp_im%val(3)%point - comp_im%val(2)%point)/(comp_im%x(3) - comp_im%x(2)) * (Xc-comp_im%x(2)) + comp_im%val(2)%point
				res = (f_up-f_low)/(comp_im%y(2)-comp_im%y(1))*(Yc-comp_im%y(1)) + f_low
			else
				if(comp_im%kas_paremvasak) then
					if(Xc<comp_im%x(5)) then
						res = get_comp_im_val(comp_im%sub1, Xc, Yc)
					else
						res = get_comp_im_val(comp_im%sub2, Xc, Yc)
					end if
				else
					if(Yc<comp_im%y(5)) then
						res = get_comp_im_val(comp_im%sub1, Xc, Yc)
					else
						res = get_comp_im_val(comp_im%sub2, Xc, Yc)
					end if
				end if
			end if
				
		end function get_comp_im_val !kontrollitud
		subroutine fill_comp_image_real(los_val_func, image_number, new) !default on see, et ei tehta uut, vaid voetakse image_number
			implicit none
			logical, intent(in), optional :: new
			integer, intent(inout) :: image_number
			real(rk), dimension(:), pointer :: val_ladu
			real(rk), dimension(:), pointer :: x_ladu
			real(rk), dimension(:), pointer :: y_ladu
			logical, dimension(:), pointer :: kas_arvutatud_ladu
			integer, dimension(:), pointer :: id_list
			type(comp_image_type), pointer :: juurikas !puu alus, kust otsimisi peaks tegema hakkama
			integer :: ladu_jrk !loendab, ladu indekseid
			type(comp_image_linked_list_type), pointer :: suur_pilt, uusjupp
			integer :: i
			interface
				function los_val_func(Xc, Yc) result(res)
					import rk
					real(rk) :: res
					real(rk), intent(in) :: Xc, Yc
				end function los_val_func
			end interface
			
			!
			! ========= oige comp_im otsimine linked listist =========
			!
			suur_pilt=>comp_im_list
			if(present(new) .and. .not.new) then !ehk kui pole vaja uut teha
				do i=1,im_list_maxsize
					if(suur_pilt%mitmes_see_comp_im_on == image_number) then
						exit !ehk see on oige pilt
					end if
					if(.not.associated(suur_pilt%next)) then
						stop("otsitud comp_im pilti pole olemas")
					end if
					suur_pilt=>suur_pilt%next
				end do
			else !lisatakse uus linked listi	
				do while(associated(suur_pilt%next))
					suur_pilt=>suur_pilt%next
				end do
				allocate(uusjupp); 
				uusjupp%mitmes_see_comp_im_on = suur_pilt%mitmes_see_comp_im_on + 1
				suur_pilt%next => uusjupp
				image_number = suur_pilt%mitmes_see_comp_im_on
			end if
			!
			! ============ algne tyybi elementide initsialiseerimine ==========
			!
			juurikas => suur_pilt%comp_im !p2ris see pilt
			ladu_jrk = 1 !algv22rtus loendajale
			allocate(x_ladu(maxsize)); x_ladu=>suur_pilt%x_ladu
			allocate(y_ladu(maxsize)); y_ladu=>suur_pilt%y_ladu
			allocate(id_list(maxsize)); id_list=>suur_pilt%id_list
			allocate(kas_arvutatud_ladu(maxsize)); kas_arvutatud_ladu=>suur_pilt%kas_arvutatud_ladu
			allocate(val_ladu(1:maxsize)); val_ladu=>suur_pilt%val_ladu
			x_ladu = 0
			y_ladu = 0
			val_ladu(:) = 0
			kas_arvutatud_ladu = .false.
	
! 			======= reaalne arvutamine ======== 
			
			call fill_im(suur_pilt%comp_im, x0_default, y0_default, x1_default, y1_default, 0)


		contains
			recursive subroutine fill_im(res, x0, y0, x1, y1, mis_levelil)
				implicit none
				real(rk), intent(in) :: x0, y0, x1, y1
				type(comp_image_type), intent(inout) :: res
				real(rk) :: xc, yc !abimuutuja... kasti tsenter
				real(rk) :: pindala !abimuutuja
				integer, intent(in) :: mis_levelil !testmuutuja ma max leveliga vordlemiseks
				logical :: kas_juba_olemas !tulemuse kontroll, kas varem juba arvutatud
				integer :: id_juba_olemas
				integer :: i
				integer, save :: vidin_counter = 0
				
				! ==================== gridi raku t2itmine ==================== 
				xc = 0.5*(x0+x1); yc = 0.5*(y0+y1)
				res%x = [x0, x0, x1, x1, xc];  
				res%y = [y0, y1, y1, y0, yc]
				res%id = puudub
				do i=1,5
					vidin_counter = vidin_counter + 1
					if(.not.associated(res%val(i)%point)) then
! 						call kas_varem_arvutatud0(res%x(i), res%y(i), kas_juba_olemas, id_juba_olemas) !kohutavalt aegland, aga toimib
! 						call kas_varem_arvutatud1(res%x(i), res%y(i), kas_juba_olemas, id_juba_olemas) !alati uuesti arvutada
						call kas_varem_arvutatud(res%x(i), res%y(i), kas_juba_olemas, id_juba_olemas)
						if(kas_juba_olemas) then
							res%val(i)%point => val_ladu(id_juba_olemas)
							res%id(i) = id_juba_olemas
						else
							call lisa_xy_lattu(res%x(i), res%y(i))
							res%val(i)%point => val_ladu(ladu_jrk-1)
							res%id(i) = ladu_jrk-1
						end if						
					end if
				end do
				
!
! 				! ==================== alamstruktuuri vajadus ====================
!
				res%last_level = .not.kas_jagada_edasi(res%x, res%y, res%val)  !tavaline kriteerium
				res%last_level = res%last_level .or. (mis_levelil>comp_im_maxlevel)  !teatud sygavustesse ei tohi minna
				res%last_level = res%last_level .and. .not.(mis_levelil < 2) !v2hemalt 2 levelit peab olema
!
! 				! ==================== alamstruktuuri t2itmine ====================
!				
				if(associated(res%sub1) .or. associated(res%sub2)) stop("probleemid sub1 voi sub2 associated-ga")
				if(.not.res%last_level) then !kui on subgrid t2psust vaja...
					res%kas_paremvasak = leia_jagamise_suund(res%x, res%y, res%val) !jagatakse vastavalt sellele, mis suunas on tohusam...
					if(res%kas_paremvasak) then
						!arvutatakse rekrusiivselt edasi parem ja vasak pool
						allocate(res%sub1)
						call fill_im(res=res%sub1, x0=x0, y0=y0, x1=xc, y1=y1, mis_levelil=mis_levelil+1) !vasak
						 allocate(res%sub2)
						call fill_im(res=res%sub2, x0=xc, y0=y0, x1=x1, y1=y1, mis_levelil=mis_levelil+1) !parem
					else
						!alumine on enne ylemist... nii, et koordinaat alati kasvab
						allocate(res%sub1)
						call fill_im(res=res%sub1, x0=x0, y0=y0, x1=x1, y1=yc, mis_levelil=mis_levelil+1)  !alumine
					 	allocate(res%sub2)
						call fill_im(res=res%sub2, x0=x0, y0=yc, x1=x1, y1=y1, mis_levelil=mis_levelil+1) !ylemine
					end if
				end if
!
! 				! ==================== summaarsete v22rtuste t2itmine ====================
! !
! 				if(res%last_level) then
! 					pindala = (res%x(3)-res%x(1)) * (res%y(3) - res%y(1))
! 					res%sum_val%rk = 0.5*res%val(5)%point + 0.125*(res%val(1)%point+res%val(2)%point+res%val(3)%point+res%val(4)%point) !keskmine v22rtus... piksli nurgad l2hevd 1/4 kaaluga arvesse
! 					res%sum_val_korda_pind%rk = res%sum_val%rk * pindala
! 				else
! 					res%sum_val_korda_pind%rk = res%sub1%sum_val_korda_pind%rk + res%sub2%sum_val_korda_pind%rk
! 				end if
			end subroutine fill_im	
			function kas_jagada_edasi(x,y,val) result(res)
				!kontorllib, kas on piksel piisavalt sygavale arvutatud, et lineaarne interpoleerimine on sobilik.
				implicit none
				real(rk), dimension(1:5), intent(in) :: x,y
				type(rk_point_type), dimension(1:5), intent(in) :: val 
				logical :: res
				real(rk) :: val0
				
				!diagonaale pidi interpoleerib keskkohta
				val0 = abs(max(  abs(0.5*(val(1)%point+val(3)%point)-val(5)%point),   abs(0.5*(val(2)%point+val(4)%point)-val(5)%point)  ))
				res = val0 > (comp_im_edasijagamise_threshold * val(5)%point)
				res = res .and. (abs(x(3)-x(1)) > 0.001) .and. (abs(y(3)-y(1)) > 0.001) !1 pc miinimum gridi tihedus
! 				res = res .and. val0 > comp_im_edasijagamise_maksimaalne_abs_t2psus !absoluutne tiheduse piir, et v2ltida liiga pisisust
			end function kas_jagada_edasi
			subroutine kas_varem_arvutatud0(x,y,kas_varem, id_varem)	!(kohutavalt) aeglane ja lihtne testversioon
				implicit none
				real(rk), intent(in) :: x,y !see koht, mida kontrollitakse
				logical, intent(out) :: kas_varem !kontroll, kas varem arvutatud
				integer, intent(out) :: id_varem !v2ljastatakse siis kui varem arvutatud
				logical, dimension(1:maxsize) :: mask
				id_varem = -1 !ehk puudulik
				mask = kas_arvutatud_ladu .and. abs(x-x_ladu)<epsilon(0.0) .and. abs(y-y_ladu)<epsilon(0.0)
				kas_varem = any(mask)
				if(kas_varem) id_varem = maxloc(x_ladu, mask=mask, dim=1)
! kas_varem = .false.
			end subroutine kas_varem_arvutatud0
			subroutine kas_varem_arvutatud1(x,y,kas_varem, id_varem)	!alati arvutab uuesti
				implicit none
				real(rk), intent(in) :: x,y !see koht, mida kontrollitakse
				logical, intent(out) :: kas_varem !kontroll, kas varem arvutatud
				integer, intent(out) :: id_varem !v2ljastatakse siis kui varem arvutatud
				logical, dimension(1:maxsize) :: mask
				kas_varem = .false.
				id_varem = puudub
! kas_varem = .false.
			end subroutine kas_varem_arvutatud1
			subroutine kas_varem_arvutatud(x,y,kas_varem, id_varem)		!otsib alati mooda puud oige v22rtuse
				implicit none
				real(rk), intent(in) :: x,y !see koht, mida kontrollitakse
				logical, intent(out) :: kas_varem !kontroll, kas varem arvutatud
				integer, intent(out) :: id_varem 
				kas_varem = .false.
				call otsi(juurikas, x,y, kas_varem, id_varem)
			end subroutine kas_varem_arvutatud
			recursive subroutine otsi(im, x,y, kas_varem, id_varem)		!kuulub kokku funtktsiooniga kas_varem_arvutatud
				implicit none
				real(rk), intent(in) :: x,y !see koht, mida kontrollitakse
				logical, intent(inout) :: kas_varem !kontroll, kas varem arvutatud
				integer, intent(out) :: id_varem 
				type(comp_image_type), intent(in) :: im
				logical, dimension(1:5) :: mask
				integer :: i
				
				id_varem = puudub
! 				mask = (abs(im%x-x)<epsilon(1.0)) .and. (abs(im%y-y)<epsilon(1.0))
				mask = (im%y==y .and. im%x==x)
				if(any(mask)) then
					id_varem = sum(im%id, mask = mask, dim=1) !sisuliselt ainsa v22rtuse saamine
					kas_varem = id_varem /= puudub !see on kontroll, kas kontrollib iseenda node-i
				else
					if(im%kas_paremvasak) then
						!alati tuleb sub1 enne kontrollida kui sub2 sest t2itmine toimub selles j2rjekorras
! 						if(x<im%x(5) .or. abs(im%x(5)-x)<epsilon(1.0) ) then
						if(x<im%x(5) .or. im%x(5)==x) then
							if(associated(im%sub1)) call otsi(im%sub1, x,y,kas_varem, id_varem) 
						else
							if(associated(im%sub2)) call otsi(im%sub2, x,y,kas_varem, id_varem)
						end if
					else
! 						if(y<im%y(5) .or. abs(im%y(5)-y)<epsilon(1.0) ) then
						if(y<im%y(5) .or. im%y(5)==y) then
							if(associated(im%sub1)) call otsi(im%sub1, x,y,kas_varem, id_varem) !alumine
						else
							if(associated(im%sub2)) call otsi(im%sub2, x,y,kas_varem, id_varem)
						end if
					end if
				end if
			end subroutine otsi
			subroutine lisa_xy_lattu(x,y)
				implicit none
				real(rk), intent(in) :: x,y
				val_ladu(ladu_jrk) =los_val_func(x, y)
				x_ladu(ladu_jrk) = x
				y_ladu(ladu_jrk) = y
				kas_arvutatud_ladu(ladu_jrk) = .true.
				id_list(ladu_jrk) = ladu_jrk
				ladu_jrk = ladu_jrk + 1
				if(ladu_jrk>maxsize) then
					stop "comp_image: fill_im ... counter liiga suur massiivi jaoks"
				end if
! 				print*, ladu_jrk
			end subroutine lisa_xy_lattu
			function leia_jagamise_suund(x,y,val) result(res)
				implicit none
				real(rk), intent(in), dimension(1:5) :: x,y
				type(rk_point_type), dimension(5), intent(in) :: val 
				real(rk), pointer :: selecttest
				logical :: res
				real(rk), dimension(1:4) :: test, test2
				integer :: max1, max2 !ei kasuta
				logical :: kas_kasti_suunas, tmpl

					!kontrollib, kus suunas standardh2lve suurim ning valib selle... vordse korral poolitab pikimat pidi
					test = [val(1)%point, val(2)%point, val(3)%point, val(4)%point]
					test2 = [ (0.5*(test(1)+test(2)) ), (0.5*(test(2)+test(3)) ), (0.5*(test(3)+test(4)) ), (0.5*(test(1)+test(4)) ) ]
					kas_kasti_suunas = .not.(count(test2>val(5)%point)==2) !ehk kas tegelik v22rtus on nurkade v22rtuste keskel (suuruselt 2. ja 3. vahel)... kui ei ole, siis teeb ruudu sarnasemaks
					
					!j2rgnevad eraldavad suurima ja suuruselt teise hajumisega suuna
					test = abs(test-val(5)%point)
					max1 = maxloc(test, dim=1); test(max1) = -1.0; max2 = maxloc(test, dim=1)
					!kontrollib, kas suurima muutuse suund on vastaskylgedes
! 					if(	kas_kasti_suunas .or. (max1==1 .and. max2==3) .or. (max1==3 .and. max2==1) .or. (max1==2 .and. max2==4) .or. (max1==4 .and. max2==2)) then
					if(	kas_kasti_suunas .or. (abs(max1 - max2)==2)) then 
						res = abs(x(1)-x(3)) > abs(y(1)-y(3))
					else
						! res = (max1==1 .and. max2==2) .or. (max1==2 .and. max2==1) .or. (max1==3 .and. max2==4) .or. (max1==4 .and. max2==3)
						! res = .not.res
						res = (max1<3 .and. max2<3) .or. (max1>2 .and. max2>2)
						res = .not.res
					end if
					
			end function leia_jagamise_suund
		end subroutine fill_comp_image_real
		subroutine kirjuta_comp_im_faili(comp_im_number, file_name)
			implicit none
			integer, intent(in) :: comp_im_number
			type(comp_image_type), pointer :: comp_im
			character(len=*), intent(in) :: file_name
			integer :: iunit
			call get_pointer_to_comp_im_number_X(comp_im, comp_im_number)
			open(unit=iunit, file=file_name, action="write")
			call kirjuta_kast(comp_im)
			close(iunit)
		contains
			recursive subroutine kirjuta_kast(comp_im)
				implicit none
				type(comp_image_type), intent(in) :: comp_im
				integer :: i
				do i=1,5
					write(unit=iunit, fmt = "(3F16.8)") comp_im%x(i), comp_im%y(i), comp_im%val(i)%point
				end do
				if(.not.comp_im%last_level) then
					call kirjuta_kast(comp_im%sub1)
					call kirjuta_kast(comp_im%sub2)
				end if
			end subroutine kirjuta_kast
		end subroutine kirjuta_comp_im_faili
		subroutine get_pointer_to_comp_im_number_X(comp_im, number)
			implicit none
			type(comp_image_type), pointer, intent(out) :: comp_im
			integer, intent(in) :: number
			integer :: i
			type(comp_image_linked_list_type), pointer :: cil
! 							print*, "tahab im number", number
			cil => comp_im_list
			do i=1,im_list_maxsize
				if(cil%mitmes_see_comp_im_on == number) then
					exit !ehk see on oige pilt
				end if
				if(.not.associated(cil%next)) then
					stop("otsitud comp_im pilti pole olemas ....")
				end if
				cil=>cil%next
			end do
			comp_im => cil%comp_im
		end subroutine get_pointer_to_comp_im_number_X
		
end module comp_image_real_module