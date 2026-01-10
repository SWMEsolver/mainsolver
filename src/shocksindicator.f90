!===============================================================================!
MODULE MOD_ShocksIndicator
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ShocksIndicatorX
  MODULE PROCEDURE ShocksIndicatorX
END INTERFACE

INTERFACE ShocksIndicatorY
  MODULE PROCEDURE ShocksIndicatorY
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: ShocksIndicatorX
PUBLIC :: ShocksIndicatorY
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!


!===============================================================================!
SUBROUTINE ShocksIndicatorX()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: V
USE MOD_FiniteVolume_vars,ONLY: nVar, nXs, nYs, Ind
IMPLICIT NONE
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj
LOGICAL :: Shock

Shock      = .FALSE.
Ind(1,:,:) = .FALSE.

DO jj = 1,nYs
    DO ii = 1,nXs
        Shock = JamesonIndicator(V(1:nVar,ii-1:ii+1,jj))
        Ind(1,ii-1:ii+1,jj) = Shock
    END DO
END DO

END SUBROUTINE ShocksIndicatorX
!===============================================================================!


!===============================================================================!
SUBROUTINE ShocksIndicatorY()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: V
USE MOD_FiniteVolume_vars,ONLY: nXs, nYs, Ind
IMPLICIT NONE
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj
LOGICAL :: Shock
!-------------------------------------------------------------------------------!

Shock      = .FALSE.
Ind(2,:,:) = .FALSE.

DO jj = 1,nYs
    DO ii = 1,nXs
        Shock = JamesonIndicator(V(:,ii,jj-1:jj+1))
        Ind(2,ii,jj-1:jj+1) = Shock
    END DO
END DO

END SUBROUTINE ShocksIndicatorY
!===============================================================================!


!===============================================================================!
FUNCTION JamesonIndicator(Prim)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: Prim(1:nVar,-1:1)
LOGICAL         :: JamesonIndicator
!-------------------------------------------------------------------------------!
REAL            :: Indicator
REAL            :: IndicatorLimit
REAL            :: IndVar(-1:1)


IndicatorLimit = 0.5E-03

JamesonIndicator = .FALSE.

IndVar    = Prim(1,:)
Indicator = ABS(IndVar(-1) + IndVar(+1) - 2.0*IndVar(0))
Indicator = Indicator/(ABS(IndVar(-1)) + ABS(IndVar(+1)) + ABS(2.0*IndVar(0)))

IF (Indicator .GT. IndicatorLimit) THEN
    JamesonIndicator = .TRUE.
END IF

END FUNCTION JamesonIndicator
!===============================================================================!


!-------------------------------------------------------------------------------!
END MODULE MOD_ShocksIndicator
!===============================================================================!
