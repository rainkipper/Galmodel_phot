!*********************
! @author elmo tempel
! @date 22.03.2011
! @lines 13
!*********************
! defines real type - normal or double
Module realkind
  Implicit None
  integer,parameter  :: rk = kind(1.0D0)  !RealKind: double
  !integer,parameter  :: rk = kind(1.0)    !normal precision
  integer,parameter  :: rsp = kind(1.0)  !RealKind: single
  integer,parameter  :: rdp = kind(1.0D0)  !RealKind: double
End Module realkind
