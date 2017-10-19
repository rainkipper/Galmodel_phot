!*********************
! @author elmo tempel
! @date 14.10.2011
! @lines 129
!*********************
! see on koopia nr amoeba moodulist koos koigi lisadega
! eraldiseisva moodulina saan kasutada oma rk kind parameetrit
! ning kiiruses on ehk ka mingi efekt
MODULE amoeba_module
  use constants_module
  implicit none
  private
  public:: my_amoeba
contains
  ! amoeba fitting subroutine, taken from NR
  SUBROUTINE my_amoeba(p,y,ftol,func,iter)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: iter !n2itab, mitu iteratsiooni tegi
    REAL(rk), INTENT(IN) :: ftol !fractional convergence tolerance
    REAL(rk), DIMENSION(:), INTENT(INOUT) :: y ! (N+1) punktile likelihoodi vm v22rtused, mida minimeerida.. ehk likelihood tuleb teisipidi defineerida
    REAL(rk), DIMENSION(:,:), INTENT(INOUT) :: p !(N+1,N) matrix, esimene on puntid, teine parameetrite v22rtused
    INTERFACE
       FUNCTION func(x)
         import rk
         IMPLICIT NONE
         REAL(rk), DIMENSION(:), INTENT(IN) :: x
         REAL(rk) :: func
       END FUNCTION func
    END INTERFACE
    INTEGER, PARAMETER :: ITMAX=5000
    REAL(rk), PARAMETER :: TINY=1.0e-12
    INTEGER :: ihi,ndim
    REAL(rk), DIMENSION(SIZE(p,2)) :: psum
    CALL amoeba_private
  CONTAINS
    FUNCTION imaxloc(arr)
      REAL(rk), DIMENSION(:), INTENT(IN) :: arr
      INTEGER, DIMENSION(1) :: imax
      INTEGER :: imaxloc
      imax=MAXLOC(arr(:))
      imaxloc=imax(1)
    END FUNCTION imaxloc
    FUNCTION iminloc(arr)
      REAL(rk), DIMENSION(:), INTENT(IN) :: arr
      INTEGER, DIMENSION(1) :: imin
      INTEGER :: iminloc
      imin=MINLOC(arr(:))
      iminloc=imin(1)
    END FUNCTION iminloc
    SUBROUTINE swap_r(a,b)
      REAL(rk), INTENT(INOUT) :: a,b
      REAL(rk) :: dum
      dum=a
      a=b
      b=dum
    END SUBROUTINE swap_r
    SUBROUTINE swap_rv(a,b)
      REAL(rk), DIMENSION(:), INTENT(INOUT) :: a,b
      REAL(rk), DIMENSION(SIZE(a)) :: dum
      dum=a
      a=b
      b=dum
    END SUBROUTINE swap_rv
    SUBROUTINE amoeba_private
      IMPLICIT NONE
      INTEGER :: i,ilo,inhi
      REAL(rk) :: rtol,ysave,ytry,ytmp
      ndim=SIZE(p,2)
      iter=0
      psum(:)=SUM(p(:,:),dim=1)
      DO
         ilo=iminloc(y(:))
         ihi=imaxloc(y(:))
         ytmp=y(ihi)
         y(ihi)=y(ilo)
         inhi=imaxloc(y(:))
         y(ihi)=ytmp
         rtol=2.0_rk*ABS(y(ihi)-y(ilo))/(ABS(y(ihi))+ABS(y(ilo))+TINY)
         IF (rtol < ftol) THEN
            CALL swap_r(y(1),y(ilo))
            CALL swap_rv(p(1,:),p(ilo,:))
            RETURN
         END IF
         IF (iter >= ITMAX) THEN
            PRINT*, "ITMAX exceeded in amoeba..."
            CALL swap_r(y(1),y(ilo))
            CALL swap_rv(p(1,:),p(ilo,:))
            RETURN
            print*, "siin ei tohi olla... itmax ex amoeba..."
         END IF
         ytry=amotry(-1.0_rk)
         iter=iter+1
         IF (ytry <= y(ilo)) THEN
            ytry=amotry(2.0_rk)
            iter=iter+1
         ELSE IF (ytry >= y(inhi)) THEN
            ysave=y(ihi)
            ytry=amotry(0.5_rk)
            iter=iter+1
            IF (ytry >= ysave) THEN
               p(:,:)=0.5_rk*(p(:,:)+SPREAD(p(ilo,:),1,SIZE(p,1)))
               DO i=1,ndim+1
                  IF (i /= ilo) y(i)=func(p(i,:))
               END DO
               iter=iter+ndim
               psum(:)=SUM(p(:,:),dim=1)
            END IF
         END IF
      END DO
    END SUBROUTINE amoeba_private
    FUNCTION amotry(fac)
      IMPLICIT NONE
      REAL(rk), INTENT(IN) :: fac
      REAL(rk) :: amotry
      REAL(rk) :: fac1,fac2,ytry
      REAL(rk), DIMENSION(SIZE(p,2)) :: ptry
      fac1=(1.0_rk-fac)/ndim
      fac2=fac1-fac
      ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
      ytry=func(ptry)
      IF (ytry < y(ihi)) THEN
         y(ihi)=ytry
         psum(:)=psum(:)-p(ihi,:)+ptry(:)
         p(ihi,:)=ptry(:)
      END IF
      amotry=ytry
    END FUNCTION amotry
  END SUBROUTINE my_amoeba
END MODULE amoeba_module
