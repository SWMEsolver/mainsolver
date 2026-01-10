!===============================================================================!
MODULE NonlinearSolver
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE NewtonRaphsonSolver
    MODULE PROCEDURE NewtonRaphsonSolver
END INTERFACE

!-------------------------------------------------------------------------------!
PUBLIC :: NewtonRaphsonSolver
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!


!===============================================================================!
SUBROUTINE NewtonRaphsonSolver(Cons0, Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)  :: Cons0(1:nVar)
REAL, INTENT(OUT) :: Cons( 1:nVar)
!-------------------------------------------------------------------------------!
REAL    :: Res(1:nVar), dCons(1:nVar), matrixJ(1:nVar,1:nVar), NormRes, Norm0
INTEGER :: iter, maxIter, info
INTEGER, ALLOCATABLE :: ipiv(:)
REAL, PARAMETER      :: tol = 1.0E-6 ! Newton tolerance
  
Cons = Cons0   ! Initial guess
maxIter = 100

DO iter = 1, maxIter

    CALL NonlinearFunction(Cons0, Cons, Res)
    Call ComputeJacobianCD(Cons0, Cons, matrixJ)
    

    NormRes = SQRT(SUM(Res*Res))
     
    IF (iter == 1) THEN
        Norm0 = MAX(NormRes, 1.0)
    END IF
    IF (NormRes < tol*Norm0) THEN
        EXIT
    END IF
    !-----------------------------------------------------------------!
    ! Solve J * dCons = -Res
    dCons = -Res
    ALLOCATE(ipiv(nVar))
    CALL dgesv(nVar, 1, matrixJ, nVar, ipiv, dCons, nVar, info)
    DEALLOCATE(ipiv)
    !-----------------------------------------------------------------!
    IF (info /= 0) THEN
        PRINT *, 'Linear solver failed, info=', info
        STOP
    END IF
    Cons = Cons + dCons
    
END DO


IF (iter > maxIter) THEN
    PRINT *, 'Newton did not converge within maxIter'
    STOP
END IF
    
END SUBROUTINE NewtonRaphsonSolver
!===============================================================================!


!===============================================================================!  
SUBROUTINE ComputeJacobianFD(Cons0, Cons, matrixJ)
USE MOD_FiniteVolume_vars,ONLY: nVar
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)  :: Cons(1:nVar), Cons0(1:nVar)
REAL, INTENT(OUT) :: matrixJ(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
REAL    :: SS0(1:nVar), SSplus(1:nVar), ConsTemp(1:nVar), EPSFD
INTEGER :: iin, jjn

SS0   = 0.0; SSplus = 0.0; ConsTemp = 0.0
EPSFD = 1.0E-8 !tolerance
ConsTemp = Cons

CALL NonlinearFunction(Cons0, Cons, SS0)

DO jjn = 1, nVar

    ConsTemp      = Cons
    ConsTemp(jjn) = Cons(jjn) + EPSFD
    CALL NonlinearFunction(Cons0, ConsTemp, SSplus)
    
    DO iin = 1, nVar
        matrixJ(iin, jjn) = (SSplus(iin) - SS0(iin))/EPSFD
    END DO
END DO
END SUBROUTINE ComputeJacobianFD
!===============================================================================!


!===============================================================================!
SUBROUTINE ComputeJacobianCD(Cons0, Cons, matrixJ)
USE MOD_FiniteVolume_vars,ONLY: nVar
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)  :: Cons(1:nVar), Cons0(1:nVar)
REAL, INTENT(OUT) :: matrixJ(1:nVar,1:nVar)
!-------------------------------------------------------------------------------!
REAL    :: SSplus(1:nVar), SSminus(1:nVar), ConsTemp(1:nVar), EPSCD
INTEGER :: iin, jjn

SSminus = 0.0; SSplus = 0.0; ConsTemp = 0.0
EPSCD   = 1.0E-8 ! tolerance

DO jjn = 1, nVar

    ConsTemp      = Cons
    ConsTemp(jjn) = Cons(jjn) + EPSCD
    CALL NonlinearFunction(Cons0, ConsTemp, SSplus)
    
    ConsTemp      = Cons
    ConsTemp(jjn) = Cons(jjn) - EPSCD
    CALL NonlinearFunction(Cons0, ConsTemp, SSminus)

    DO iin = 1, nVar
        matrixJ(iin, jjn) = (SSplus(iin) - SSminus(iin))/(2.0*EPSCD)
    END DO
END DO

END SUBROUTINE ComputeJacobianCD
!===============================================================================!


!===============================================================================!
SUBROUTINE NonlinearFunction(Cons0, Cons, SSout)
USE MOD_FiniteVolume_vars,ONLY: nVar, dt
USE MOD_MomentModels,     ONLY: SourceSplit
USE MOD_PhysicsFrame,     ONLY: ConsToPrim
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)  :: Cons(1:nVar), Cons0(1:nVar)
REAL, INTENT(OUT) :: SSout(1:nVar)
!-------------------------------------------------------------------------------!
REAL    :: Prim(1:nVar)
SSout = 0.0
!-------------------------------------------------------------------------------!
! du/dt = S(u)  ==> u^{n+1} = u^n + dt*S(u^{n+1})
!               ==> (u^{n+1} - dt*S(u^{n+1})) - u^n = 0
!-------------------------------------------------------------------------------!
CALL ConsToPrim(Cons, Prim)
CALL SourceSplit(Prim, SSout)
SSout = Cons - Cons0 - dt*SSout 
!-------------------------------------------------------------------------------!
END SUBROUTINE NonlinearFunction
!===============================================================================!

END MODULE NonlinearSolver
!===============================================================================!
