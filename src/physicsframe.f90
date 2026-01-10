!===============================================================================!
MODULE MOD_PhysicsFrame
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ExactFunction
    MODULE PROCEDURE ExactFunction
END INTERFACE

INTERFACE ExactFunctionWB
    MODULE PROCEDURE ExactFunctionWB
END INTERFACE

INTERFACE BoundaryConditions
    MODULE PROCEDURE BoundaryConditions
END INTERFACE

INTERFACE ConsToPrim
    MODULE PROCEDURE ConsToPrim
END INTERFACE

INTERFACE PrimToCons
    MODULE PROCEDURE PrimToCons
END INTERFACE

INTERFACE Bathymetry   
    MODULE PROCEDURE Bathymetry   
END INTERFACE

INTERFACE DerivativeBathymetry
    MODULE PROCEDURE DerivativeBathymetry
END INTERFACE

INTERFACE Bathymetry_X   
    MODULE PROCEDURE Bathymetry_X
END INTERFACE

INTERFACE Bathymetry_Y
    MODULE PROCEDURE Bathymetry_Y
END INTERFACE

INTERFACE RootsLegendre01
    MODULE PROCEDURE RootsLegendre01
END INTERFACE

INTERFACE RootsLegendre
    MODULE PROCEDURE RootsLegendre
END INTERFACE

INTERFACE LegendrePoly
    MODULE PROCEDURE LegendrePoly
END INTERFACE

INTERFACE LegendrePolyNormalized01
    MODULE PROCEDURE LegendrePolyNormalized01
END INTERFACE

INTERFACE JacobiRoots
    MODULE PROCEDURE JacobiRoots
END INTERFACE

!-------------------------------------------------------------------------------!
PUBLIC :: ExactFunction
PUBLIC :: ExactFunctionWB
PUBLIC :: BoundaryConditions
PUBLIC :: ConsToPrim
PUBLIC :: PrimToCons
PUBLIC :: Bathymetry
PUBLIC :: DerivativeBathymetry
PUBLIC :: Bathymetry_X
PUBLIC :: Bathymetry_Y

PUBLIC :: RootsLegendre01
PUBLIC :: RootsLegendre
PUBLIC :: LegendrePoly
PUBLIC :: LegendrePolyNormalized01

PUBLIC :: JacobiRoots
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!


!===============================================================================!
SUBROUTINE ExactFunction(t,x,Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars, ONLY: nVar, nDims, nMoms
USE MOD_FiniteVolume_vars, ONLY: MESH_X0, MESH_SX, Hmin
USE MOD_FiniteVolume_vars, ONLY: PI, Gravity
USE MOD_FiniteVolume_vars, ONLY: InitialFlag, VarInitOpt
USE MOD_FiniteVolume_vars, ONLY: ParameterLambda, ParameterNu
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: t
REAL,INTENT(IN)    :: x(1:nDims)
REAL,INTENT(OUT)   :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
REAL               :: Prim(1:nVar)
REAL               :: xc(1:nDims), xm(1:nDims), r, r0, r2, aux_vinf(1:2)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!
!*OUR VARIABLES
REAL               :: nu, lmbda
REAL               :: g
nu = ParameterNu
lmbda = ParameterLambda
g = Gravity
Cons = 0.0
Prim = 0.0

xc = 0.0
xm = 0.0 + 0.0*t

aux_vinf = 0.0
SELECT CASE (InitialFlag)
CASE (1) 
    VarInitOpt = "Exponentially curved h(x)"
    ! Case JK 2020:
    Prim(1) =  1.0 + exp(3.0*cos(PI*(x(1) + 0.5)))/exp(4.0)      

    ! Reference simulation, g = 1.0, lambda = nu = 0.1, t = 2
    ! Linear case u = 0.5*zeta
    Prim(2) =  0.25 ! um
    Prim(3) = -0.25 ! alpha1
    ! Quadratic case, u = 1.5*zeta - 1.5*zeta**2
    !Prim(2) =  0.25 ! um
    !Prim(3) =  0    ! alpha1
    !Prim(4) = -0.25 ! alpha2

CASE (2) !*LAKE AT REST
    VarInitOpt = "Lake-at-rest"
    Prim(1) = 1. - 0.1*SIN(2.*PI*x(1))*COS(2.*PI*x(nDims))
    Prim(2:nDims+1) = 0.

CASE (3) ! *PERTURBATION ANALYSIS ON WET-DRY LAKE AT REST
    VarInitOpt = "Perturbation analysis on wet-dry"
    Prim(1) = MAX(0.7 - Bathymetry(x), Hmin)
    r2 = ((x(1) + 2)**2 + (REAL(nDims) - 1.0)*(x(nDims) - 0.5)**2)*9.

    IF (r2 < 1) THEN
        Prim(1) = Prim(1) + 0.05*EXP(1.-1./(1 - r2)**2.)
    END IF
    Prim(2:nDims+1) = 0.

CASE (4) ! *CIRCULAR DAM BREAK 1
    VarInitOpt = "Circular-Dam-Break 1"
    xm = MESH_X0(1:nDims) + 0.5*MESH_SX(1:nDims)
    xc(1:nDims) = x(1:nDims) - xm(1:nDims)
    !----------------------------------------------------
    r = SQRT(xc(1)**2 + (REAL(nDims) - 1.0)*xc(nDims)**2)
    !----------------------------------------------------
    r0 = 7.
    Prim(1) = Hmin
    Prim(2:nDims+1) = 0.0
    IF (r .LE. r0) THEN
        Prim(1) = 2.5
        Prim(2:nDims+1) = 0.0
    END IF

CASE (5) ! *CIRCULAR DAM BREAK 2
    VarInitOpt = "Circular-Dam-Break 2"
    xm = MESH_X0(1:nDims) + 0.5*MESH_SX(1:nDims)
    xc(1:nDims) = x(1:nDims) - xm(1:nDims)
    r = SQRT(xc(1)**2 + (REAL(nDims) - 1.0)*xc(nDims)**2)
    r0 = 15.

    Prim(1) = 0.5
    Prim(2:nDims+1) = 0.0
    IF (r .LE. r0) THEN
        Prim(1) = 10.
        Prim(2:nDims+1) = 0.0
    END IF

CASE (6) ! *WAVE OVER DRY ISLAND
    VarInitOpt = "Wave-dry-island" 
    r2 = (x(1) + 2)**2.
    Prim(1) = MAX(0.7 - Bathymetry(x), Hmin)
    IF (r2 < 1) THEN
        Prim(1) = Prim(1) + 0.5*EXP(1.-1./(1 - r2)**2.)
    END IF
    IF (Prim(1) > Hmin*10) THEN
        Prim(2) = 1.0
    ELSE
        Prim(2) = 0.0
    END IF
    Prim(nDims+1) = 0.0

CASE (7) ! *WAVE OVER DRY ISLAND
    VarInitOpt = "Huang7" 
    Prim(1) = 1.0
    IF (x(1) < 0 ) THEN
        Prim(1) = Prim(1) + 0.5
    END IF
    Prim(2) = 0.25
    IF (nMoms .GT. 0) THEN
        Prim(3) = -0.1
    END IF
    IF (nMoms .GT. 1) THEN
        Prim(4) = -0.1
    END IF
    
CASE (8) 
    VarInitOpt = "Experiments for reduced I"
    Prim(1) = 100.*(1.0 + exp(3.0*cos(PI*(x(1) + 0.5)))/exp(4.0))
    Prim(2) = 0.0

CASE (9) 
    VarInitOpt = "Experiments for reduced II"
    Prim(1) = 0.1 *(1.0 + exp(3.0*cos(2.*PI*(x(1) + 0.5)))/exp(4.0))
    Prim(2) = 0.25

CASE (10) 
    VarInitOpt = "Experiments for inclined plane-source terms"
    Prim(1) =  0.6 + 0.2*EXP( -((X(1) - 0.5)/0.1)**2 )
    Prim(2) = 0.25
   
CASE(11)
    VarInitOpt = "Experiments for inclined plane-source terms"
    
    IF ((X(1) .GT. 0.3) .AND. (X(1) .LT. 0.5))THEN   
        Prim(1) = 0.08      
    ELSE
        Prim(1) = 1.0E-6
    END IF 
    Prim(2) =  0.25 ! um
    Prim(3) = -0.25 ! alpha1  

CASE (12) 
    VarInitOpt = "Experiments for inclined plane-source terms, zero velocity"
    
    IF ((X(1) .GT. 0.3) .AND. (X(1) .LT. 0.5))THEN   
        Prim(1) = 0.08   
        !Prim(1) = 0.08*SQRT(1.0 - ((X(1) - 0.4)**2.0)/(0.1*0.1) )                  
    ELSE
        Prim(1) = 1.0E-6
    END IF 
    Prim(2) = 0.0 ! um
    Prim(3) = 0.0 ! alpha1  
    
CASE(13)
    VarInitOpt = "Experiments for inclined plane-source terms"

    IF ( (X(1) .GT. -0.5) .AND. (X(1) .LT. 0.5) ) THEN 
        Prim(1) = SQRT( 0.5*0.5 - X(1)*X(1)) 
        Prim(2) =  0.25 ! um
        Prim(3) = -0.25 ! alpha1  
        
    ELSE
        Prim(1) = 1.0E-6
        Prim(2) = 0.0
        Prim(3) = 0.0
    END IF 

CASE(14)
    VarInitOpt = "Experiments for inclined plane-source terms"

    IF ( (X(1) .GT. 0.0) .AND. (X(1) .LT. 0.5) ) THEN 
        Prim(1) = SQRT( 0.5*0.5 - X(1)*X(1)) 
        Prim(2) =  0.25 ! um
        Prim(3) = -0.25 ! alpha1  
    ELSEIF ( (X(1) .LT. 0.0) ) THEN 
        Prim(1) = 0.5*EXP( -10.0*X(1)*X(1) ) 
        Prim(2) =  0.25 ! um
        Prim(3) = -0.25 ! alpha1             
    ELSE
        Prim(1) = 1.0E-15
        Prim(2) = 0.0
        Prim(3) = 0.0
    END IF 

CASE (15) 
    VarInitOpt = "Exponentially curved h(x) - Shallow water"
    Prim(1) =  1.0 + exp(3.0*cos(PI*(x(1) + 0.5)))/exp(4.0)      
    Prim(2) =  0.25 ! um

CASE (16) 
    VarInitOpt = "Exponentially curved h(x)"
    ! Case JK 2020:
    Prim(1) =  1.0 + exp(3.0*cos(PI*(x(1) + 0.5)))/exp(4.0)      
    Prim(2) =  0.25 !
    Prim(4) = -0.25 !
CASE DEFAULT
    ErrorMessage = "Exact function not specified"
    WRITE (*, *) ErrorMessage
    STOP
END SELECT


CALL PrimToCons(Prim, Cons)
CONTAINS
    
    REAL FUNCTION hSmoothAuxiliary(x)
        REAL, INTENT(IN) :: x
        hSmoothAuxiliary = 1.-0.1*exp(-1./atan(1.-x)**3.)
    END FUNCTION

    REAL FUNCTION hDerivSmoothAuxiliary(x)
        REAL, INTENT(IN) :: x
        hDerivSmoothAuxiliary = 3.*0.1*exp(1./atan(x - 1.)**3.)/(atan(x - 1.)**4*((x - 1.)**2 + 1.))
    END FUNCTION

END SUBROUTINE ExactFunction
!===============================================================================!


!===============================================================================!
SUBROUTINE ExactFunctionWB(x,Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nDims, PI, Hmin, InitialFlag
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: x(1:nDims)
REAL,INTENT(OUT)   :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
REAL               :: Prim(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage

SELECT CASE (InitialFlag)
  
    CASE(2) !*LAKE AT REST
  
        Prim(1) = 1. - 0.1 * SIN(2.*PI*x(1)) * COS(2.*PI*x(nDims))
        Prim(2:nDims+1) = 0.

    CASE(3) ! *PERTURBATION ANALYSIS ON WET-DRY LAKE AT REST 
        
        Prim(1)   = MAX(0.7 - Bathymetry(x), Hmin) 
        Prim(2:nDims+1) = 0.0

    CASE DEFAULT
        ErrorMessage = "Exact WB function not specified"
        WRITE(*,*) ErrorMessage
        STOP
END SELECT

CALL PrimToCons(Prim,Cons)

END SUBROUTINE ExactFunctionWB
!===============================================================================!


!===============================================================================!
REAL FUNCTION Bathymetry(X)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nDims, PI, BathymetryFlag, AngleIncl, MESH_X1
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)  :: X(1:nDims)
!-------------------------------------------------------------------------------!
REAL  :: r2
REAL  :: posA, posB, posMid, curv

SELECT CASE (BathymetryFlag)

CASE (1) ! used for *LAKE AT REST
    Bathymetry = 0.1 * SIN(2.*PI*X(1)) * COS(2.*PI*X(nDims))
    
CASE(2) ! used for *PERTURBATION ANALYSIS ON WET-DRY LAKE AT REST and *WAVE OVER DRY ISLAND 
    r2 = X(1)**2. + (REAL(nDims) - 1.0)*X(nDims)**2.
    IF (r2 < 1) THEN
        Bathymetry = EXP(1.0 - 1.0/(1.0 - r2))
    ELSE
        Bathymetry = 0.0
    END IF

CASE(3) ! Only 1D
    Bathymetry = 0.1*EXP( -((X(1) - 0.5)/0.1)**2 )

CASE(4)
        
    posA   = 0.5;
    posB   = 0.85;
    posMid = .5*(posA + posB);    
    curv   = tan(AngleIncl)
    
    
    Bathymetry = 0.0
    
    IF ((X(1) .LE. posB) .AND. (X(1) .GE. posA) ) THEN    
        Bathymetry = (.5*curv/(posB - posA))*(X(1) - posA)*(X(1) - posA);
    ELSEIF (X(1) .GT. posB) THEN
        Bathymetry = curv*(X(1) - posMid);
    END IF
    
CASE DEFAULT
    Bathymetry = 0.
END SELECT

END FUNCTION Bathymetry
!===============================================================================!


!===============================================================================!
REAL FUNCTION DerivativeBathymetry(X)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nDims, PI, BathymetryFlag, AngleIncl, MESH_X1
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)  :: X(1:nDims)
!-------------------------------------------------------------------------------!
REAL :: posA, posB, posMid, curv

SELECT CASE (BathymetryFlag)
CASE(3) ! Only 1D
    !Bathymetry          =      0.1*EXP( -((X(1) - 0.5)/0.1)**2 )
    DerivativeBathymetry = -2.*10.0*EXP( -((X(1) - 0.5)/0.1)**2 ) * (X(1) - 0.5)
    
CASE(4)
    posA   = 0.0;
    posB   = 0.7;
    posMid = .5*(posA + posB);    
    curv   = tan(AngleIncl)!/(MESH_X1(1) - posMid);
    
    DerivativeBathymetry = 0.0    
    IF ((X(1) .LE. posB) .AND. (X(1) .GE. posA) ) THEN    
        DerivativeBathymetry = curv/(posB - posA)*(X(1) - posA);
    ELSEIF (X(1) .GT. posB) THEN
        DerivativeBathymetry = curv;
    END IF
    
    
CASE DEFAULT
    DerivativeBathymetry = 0.
END SELECT

END FUNCTION DerivativeBathymetry
!===============================================================================!



!===============================================================================!
REAL FUNCTION Bathymetry_X(X)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars, ONLY: nDims, PI, BathymetryFlag   
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL, INTENT(IN)  :: X(1:nDims)
!-------------------------------------------------------------------------------!
REAL              :: r2

SELECT CASE (BathymetryFlag)
CASE (1) ! used for *LAKE AT REST
    Bathymetry_X = 0.1 * 2. * PI * COS(2.*PI*X(1)) * COS(2.*PI*X(nDims)) 

CASE(2) ! used for *PERTURBATION ANALYSIS ON WET-DRY LAKE AT REST and *WAVE OVER DRY ISLAND
    
    r2 = X(1)**2. + (REAL(nDims) - 1.0)*X(nDims)**2.
    
    IF (r2 < 1) THEN
        Bathymetry_X = -2.0*X(1)/(1.0 - r2)**2.*EXP(1.0 - 1.0/(1.0 - r2))
    ELSE
        Bathymetry_X = 0.0
    ENDIF

CASE DEFAULT
    Bathymetry_X = 0.

END SELECT

END FUNCTION Bathymetry_X
!===============================================================================!


!===============================================================================!
REAL FUNCTION Bathymetry_Y(X)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars, ONLY: nDims, PI, nDims, BathymetryFlag   
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL, INTENT(IN) :: X(1:nDims)
!-------------------------------------------------------------------------------!
REAL             :: r2

SELECT CASE (BathymetryFlag)
CASE (1) ! used for *LAKE AT REST
    Bathymetry_Y = - 0.1 * 2. * PI * SIN(2.*PI*X(1)) * SIN(2.*PI*X(nDims))
CASE(2) ! used for *PERTURBATION ANALYSIS ON WET-DRY LAKE AT REST and *WAVE OVER DRY ISLAND
    r2 = X(1)**2. + (REAL(nDims) - 1.0)*X(nDims)**2.
    
    IF (r2 < 1) THEN
        Bathymetry_Y = -2.0*X(nDims)/(1.0 - r2)**2.*EXP(1.0 - 1.0/(1.0 - r2))
    ELSE
        Bathymetry_Y = 0.0
    ENDIF
CASE DEFAULT
    Bathymetry_Y = 0.
END SELECT

END FUNCTION Bathymetry_Y
!===============================================================================!


!===============================================================================!
SUBROUTINE BoundaryConditions(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars, ONLY: nVar, nXs, nYs, nGhostsX, nGhostsY
USE MOD_FiniteVolume_vars, ONLY: PrimRefState1, PrimRefState2 
USE MOD_FiniteVolume_vars, ONLY: PrimRefState3, PrimRefState4
USE MOD_FiniteVolume_vars, ONLY: BoundaryFlag, VarBConOpt
USE MOD_FiniteVolume_vars, ONLY: U, V
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: t
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: xt
REAL               :: Prim_in(1:nVar), Prim_out(1:nVar)
REAL               :: Cons_in(1:nVar), Cons_out(1:nVar)
!REAL               :: ConsRefState1(1:nVar), ConsRefState2(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage
xt = 0.0*t
idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!

PrimRefState1 = 0.0
PrimRefState2 = 0.0
PrimRefState3 = 0.0
PrimRefState4 = 0.0

SELECT CASE (BoundaryFlag(1))
CASE (1) ! Periodic
    VarBConOpt(1) = "Periodic"
    DO jj = 1, nYs
        DO ii = 0, nGhostsX
            U(:, -nGhostsX + ii, jj) = U(:, nXs - nGhostsX + ii, jj)
            CALL ConsToPrim(U(:, -nGhostsX + ii, jj), V(:, -nGhostsX + ii, jj)) 
        END DO
    END DO

CASE (2) ! Transmissive
    VarBConOpt(1) = "Transmissive"
    DO jj = 1, nYs
        DO ii = 0, nGhostsX
            U(:, -nGhostsX + ii, jj) = U(:, nGhostsX - ii + 1, jj)
            CALL ConsToPrim(U(:, -nGhostsX + ii, jj), V(:, -nGhostsX + ii, jj))
        END DO
    END DO

CASE (3) ! Inflow
    VarBConOpt(1) = "Inflow"
    Prim_in = PrimRefState1
    CALL PrimToCons(Prim_in, Cons_in)
    DO jj = 1, nYs
        DO ii = 0, nGhostsX
            U(:, -nGhostsX + ii, jj) = Cons_in
            V(:, -nGhostsX + ii, jj) = PrimRefState1
            !CALL ConsToPrim(U(:, -nGhostsX + ii, jj),V(:, -nGhostsX + ii, jj))
        END DO
    END DO

CASE (4) ! Outflow
    VarBConOpt(1) = "Outflow"
    Prim_out = PrimRefState1
    CALL PrimToCons(Prim_out, Cons_out)
    DO jj = 1, nYs
        DO ii = 0, nGhostsX
            U(:, -nGhostsX + ii, jj) = Cons_out
            CALL ConsToPrim(U(:, -nGhostsX + ii, jj),V(:, -nGhostsX + ii, jj))
        END DO
    END DO

CASE (5) ! Reflecting
    VarBConOpt(1) = "Reflecting"
    DO jj = 1, nYs
        DO ii = 0, nGhostsX
            U(:, -nGhostsX + ii, jj)      =  U(:, nGhostsX - ii + 1, jj)
            U(idx_vx, -nGhostsX + ii, jj) = -U(idx_vx, nGhostsX - ii + 1, jj)
            CALL ConsToPrim(U(:, -nGhostsX + ii, jj), V(:, -nGhostsX + ii, jj))
        END DO
    END DO
CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE (*, *) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE (BoundaryFlag(2))
CASE (1) ! Periodic
    VarBConOpt(2) = "Periodic"
    DO jj = 1, nYs
        DO ii = 1, nGhostsX + 1
            U(:, nXs + ii, jj) = U(:, ii, jj)
            CALL ConsToPrim(U(:, nXs + ii, jj), V(:, nXs + ii, jj))
        END DO
    END DO
CASE (2) ! Transmissive
    VarBConOpt(2) = "Transmissive"
    DO jj = 1, nYs
        DO ii = 1, nGhostsX + 1
            U(:, nXs + ii, jj) = U(:, nXs - ii + 1, jj)
            CALL ConsToPrim(U(:, nXs + ii, jj), V(:, nXs + ii, jj))
        END DO
    END DO
CASE (3) ! Inflow
    VarBConOpt(2) = "Inflow"
    Prim_in = PrimRefState2
    CALL PrimToCons(Prim_in, Cons_in)
    DO jj = 1, nYs
        DO ii = 1, nGhostsX + 1
            U(:, nXs + ii, jj) = Cons_in
            V(:, nXs + ii, jj) = PrimRefState2
            !CALL ConsToPrim(U(:, nXs + ii, jj), V(:, nXs + ii, jj))
        END DO
    END DO
CASE (4) ! Outflow
    VarBConOpt(2) = "Outflow"
    Prim_out = PrimRefState2
    CALL PrimToCons(Prim_out, Cons_out)
    DO jj = 1, nYs
        DO ii = 1, nGhostsX + 1
            U(:, nXs + ii, jj) = Cons_out
            CALL ConsToPrim(U(:, nXs + ii, jj), V(:, nXs + ii, jj))
        END DO
    END DO
CASE (5) ! Reflecting
    VarBConOpt(2) = "Reflecting"
    DO jj = 1, nYs
        DO ii = 1, nGhostsX + 1
            U(:, nXs + ii, jj) = U(:, nXs - ii + 1, jj)
            U(idx_vx, nXs + ii, jj) = -U(idx_vx, nXs - ii + 1, jj)
            CALL ConsToPrim(U(:, nXs + ii, jj), V(:, nXs + ii, jj))
        END DO
    END DO
CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE (*, *) ErrorMessage
    STOP
END SELECT


!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE (BoundaryFlag(3))
CASE (1) ! Periodic
    VarBConOpt(4) = "Periodic"
    DO ii = 1, nXs
        DO jj = 0, nGhostsY
            U(:, ii, -nGhostsY + jj) = U(:, ii, nYs - nGhostsY + jj)            
            CALL ConsToPrim(U(:, ii, -nGhostsY + jj), V(:, ii, -nGhostsY + jj))
        END DO
    END DO
CASE (2) ! Transmissive
    VarBConOpt(4) = "Transmissive"
    DO ii = 1, nXs
        DO jj = 0, nGhostsY
            U(:, ii, -nGhostsY + jj) = U(:, ii, nGhostsY - jj + 1)
            CALL ConsToPrim(U(:, ii, -nGhostsY + jj), V(:, ii, -nGhostsY + jj))
        END DO
    END DO
CASE (3) ! Inflow
    VarBConOpt(4) = "Inflow"
    Prim_in = PrimRefState3
    CALL PrimToCons(Prim_in, Cons_in)
    DO ii = 1, nXs
        DO jj = 0, nGhostsY
            U(:, ii, -nGhostsY + jj) = Cons_in
            V(:, ii, -nGhostsY + jj) = PrimRefState3
            !CALL ConsToPrim(U(:, ii, -nGhostsY + jj), V(:, ii, -nGhostsY + jj))
        END DO
    END DO
CASE (4) ! Outflow
    VarBConOpt(4) = "Outflow"
    Prim_out = PrimRefState3
    CALL PrimToCons(Prim_out, Cons_out)
    DO ii = 1, nXs
        DO jj = 0, nGhostsY
            U(:, ii, -nGhostsY + jj) = Cons_out
            CALL ConsToPrim(U(:, ii, -nGhostsY + jj), V(:, ii, -nGhostsY + jj))
        END DO
    END DO
CASE (5) ! Reflecting
    VarBConOpt(4) = "Reflecting"
    DO ii = 1, nXs
        DO jj = 0, nGhostsY
            U(:, ii, -nGhostsY + jj) = U(:, ii, nGhostsY - jj + 1)
            U(idx_vy, ii, -nGhostsY + jj) = -U(idx_vy, ii, nGhostsY - jj + 1)
            CALL ConsToPrim(U(:, ii, -nGhostsY + jj), V(:, ii, -nGhostsY + jj))
        END DO
    END DO
CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE (*, *) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE (BoundaryFlag(4))
CASE (1) ! Periodic
    VarBConOpt(3) = "Periodic"
    DO ii = 1, nXs
        DO jj = 1, nGhostsY + 1
            U(:, ii, nYs + jj) = U(:, ii, jj)
            CALL ConsToPrim(U(:, ii, nYs + jj), V(:, ii, nYs + jj))
        END DO
    END DO
CASE (2) ! Transmissive
    VarBConOpt(3) = "Transmissive"
    DO ii = 1, nXs
        DO jj = 1, nGhostsY + 1
            U(:, ii, nYs + jj) = U(:, ii, nYs - jj + 1)
            CALL ConsToPrim(U(:, ii, nYs + jj), V(:, ii, nYs + jj))
        END DO
    END DO
CASE (3) ! Inflow
    VarBConOpt(3) = "Inflow"
    Prim_in = PrimRefState4
    CALL PrimToCons(Prim_in, Cons_in)
    DO ii = 1, nXs
        DO jj = 1, nGhostsY + 1
            U(:, ii, nYs + jj) = Cons_in
            V(:, ii, nYs + jj) = PrimRefState4
            !CALL ConsToPrim(U(:, ii, nYs + jj), V(:, ii, nYs + jj))
        END DO
    END DO
CASE (4) ! Outflow
    VarBConOpt(3) = "Outflow"
    Prim_out = PrimRefState4
    CALL PrimToCons(Prim_out, Cons_out)
    DO ii = 1, nXs
        DO jj = 1, nGhostsY + 1
            U(:, ii, nYs + jj) = Cons_out
            CALL ConsToPrim(U(:, ii, nYs + jj), V(:, ii, nYs + jj))
        END DO
    END DO
CASE (5) ! Reflecting
    VarBConOpt(3) = "Reflecting"
    DO ii = 1, nXs
        DO jj = 1, nGhostsY + 1
            U(:, ii, nYs + jj)      =  U(:, ii, nYs - jj + 1)
            U(idx_vy, ii, nYs + jj) = -U(idx_vy, ii, nYs - jj + 1)            
            CALL ConsToPrim(U(:, ii, nYs + jj), V(:, ii, nYs + jj))
        END DO
    END DO
CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE (*, *) ErrorMessage
    STOP
END SELECT



!------------------------------------!
! Upper Corners Boundary Conditions  !
!------------------------------------!
SELECT CASE (BoundaryFlag(4))
CASE (1) ! Periodic
    DO jj = -nGhostsY, 0
        DO ii = 0, nGhostsX
            U(:, -nGhostsX + ii, jj) = U(:, nXs - nGhostsX + ii, jj)
            CALL ConsToPrim(U(:, -nGhostsX + ii, jj), V(:, -nGhostsX + ii, jj))
        END DO
    END DO

    DO jj = nYs + 1, nYs + 1 + nGhostsY
        DO ii = 0, nGhostsX
            U(:, -nGhostsX + ii, jj) = U(:, nXs - nGhostsX + ii, jj)
            CALL ConsToPrim(U(:, -nGhostsX + ii, jj), V(:, -nGhostsX + ii, jj))
        END DO
    END DO

END SELECT

!--------------------------------------!
! Right Corners Boundary Conditions    !
!--------------------------------------!
SELECT CASE (BoundaryFlag(2))

CASE (1) ! Periodic
    ! First ghost cells in j
    DO jj = -nGhostsY, 0
        DO ii = 1, nGhostsX + 1
            U(:, nXs + ii, jj) = U(:, ii, jj)
            CALL ConsToPrim(U(:, nXs + ii, jj), V(:, nXs + ii, jj))
        END DO
    END DO
    
    ! Last ghost cells in j
    DO jj = nYs + 1, nYs + 1 + nGhostsY
        DO ii = 1, nGhostsX + 1
            U(:, nXs + ii, jj) = U(:, ii, jj)
            CALL ConsToPrim(U(:, nXs + ii, jj), V(:, nXs + ii, jj))
        END DO
    END DO

END SELECT

!DO jj = 1, nYs
!    DO ii = 1, nXs
!        CALL ConsToPrim( U(:, ii, jj),  V(:, ii, jj) )
!    END DO
!END DO

END SUBROUTINE BoundaryConditions
!===============================================================================!


!===============================================================================!
SUBROUTINE ConsToPrim(Cons, Prim)
USE MOD_FiniteVolume_vars,ONLY: nVar, Hmin
IMPLICIT NONE
REAL,INTENT(IN)  :: Cons(1:nVar)
REAL,INTENT(OUT) :: Prim(1:nVar)
REAL             :: FactorMin

!---------------------------------------------------------------------------------------
! The approximation h = u_1 + min_h/u_1 is explained in https://arxiv.org/pdf/2110.13509
! see equation 32 therein.
!---------------------------------------------------------------------------------------
! Also see: doi.org/10.1016/j.amc.2015.08.121
! in page 11 (or 269 in the numeration of the paper), just below table 1.

Prim(1) = Cons(1)

IF (Prim(1) .LT. Hmin) THEN
    Prim(2:nVar) = 0.0*Prim(2:nVar)
ELSE
    FactorMin    = (2.0*Cons(1))/(Cons(1)*Cons(1) + max(Cons(1)*Cons(1),Hmin))
    Prim(2:nVar) = FactorMin*Cons(2:nVar)
END IF

END SUBROUTINE ConsToPrim
!===============================================================================!


!===============================================================================!
SUBROUTINE PrimToCons(Prim, Cons)
USE MOD_FiniteVolume_vars,ONLY: nVar
IMPLICIT NONE
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Cons(1:nVar)
Cons(1)      = Prim(1)
Cons(2:nVar) = Prim(1) * Prim(2:nVar)
END SUBROUTINE PrimToCons
!===============================================================================!


!===============================================================================!
SUBROUTINE RootsLegendre01(x, w)
USE MOD_FiniteVolume_vars,ONLY: nMoms
IMPLICIT NONE
REAL, INTENT(OUT) :: x(1:nMoms), w(1:nMoms)
REAL              :: xstd(1:nMoms), wstd(1:nMoms)
INTEGER :: ii

CALL RootsLegendre(xstd, wstd)

DO ii = 1, nMoms
   x(ii) = 0.5 * (xstd(ii) + 1.0) ! map [-1,1] → [0,1]
   w(ii) = 0.5 * wstd(ii)         ! scale weights
END DO

END SUBROUTINE RootsLegendre01
!===============================================================================!


!===============================================================================!
SUBROUTINE RootsLegendre(x, w)
! Original Gauss–Legendre on [-1,1]
USE MOD_FiniteVolume_vars,ONLY: nMoms
implicit none
REAL, INTENT(out) :: x(1:nMoms), w(1:nMoms)
INTEGER  :: ii, iter, maxiter, mMs
REAL     :: PolyN, PolyN1, dPolyN, dx, vala, tol, xi

tol     = 1e-14
maxiter = 50

mMs = (nMoms + 1)/2

DO ii = 1, mMs
    vala = real(4*ii - 1) / real(4*nMoms + 2)
    xi   = cos(acos(-1.0) * vala)

    iter = 0
    DO
        CALL LegendrePoly(xi, PolyN, dPolyN, nMoms)
        IF (abs(dPolyN) == 0) EXIT
        
        dx   = -PolyN/dPolyN
        xi   = xi + dx
        iter = iter + 1
            
        IF (abs(dx) < tol)   EXIT
        IF (iter >= maxiter) EXIT
    END DO

    x(ii) =  xi
    x(nMoms+1-ii) = -xi
    w(ii) = 2.0/((1.0 - xi*xi) * dPolyN*dPolyN)
    w(nMoms+1-ii) = w(ii)
END DO

IF (mod(nMoms,2) == 1) THEN
    x(mMs) = 0.0
END IF

END SUBROUTINE RootsLegendre
!===============================================================================!


!===============================================================================!
SUBROUTINE LegendrePoly(x, PolyN, dPolyN, nDegree)
USE MOD_FiniteVolume_vars,ONLY: nMoms
implicit none
REAL, INTENT(in)    :: x
INTEGER, INTENT(in) :: nDegree
REAL, INTENT(out)   :: PolyN, dPolyN
REAL    :: Polykm2, Polykm1, Polyk
REAL    :: dPolykm2, dPolykm1, dPolyk
INTEGER :: kk

IF (nDegree == 0) THEN
    PolyN  = 1.0
    dPolyN = 0.0
    RETURN
END IF

Polykm2 = 1.0
Polykm1 = x

dPolykm2 = 0.0
dPolykm1 = 1.0


IF (nDegree == 1) THEN
    PolyN  = Polykm1
    dPolyN = 1.0
    RETURN
END IF

do kk = 2, nDegree
    Polyk   = (REAL(2*kk-1)*x*Polykm1 - REAL(kk-1)*Polykm2) / real(kk)
    dPolyk  = REAL(2*kk-1)*Polykm1 + dPolykm2
    
    Polykm2 = Polykm1
    Polykm1 = Polyk
    
    dPolykm2 = dPolykm1
    dPolykm1 = dPolyk
        
end do

PolyN  = Polyk
dPolyN = dPolyk

END SUBROUTINE LegendrePoly
!===============================================================================!


!===============================================================================!
SUBROUTINE LegendrePolyNormalized01(x, PolyN, dPolyN, nDegree)
USE MOD_FiniteVolume_vars,ONLY: nMoms
implicit none
REAL, INTENT(in)    :: x
INTEGER, INTENT(in) :: nDegree
REAL, INTENT(out)   :: PolyN, dPolyN
REAL    :: Polykm2, Polykm1, Polyk
INTEGER :: kk

CALL LegendrePoly(2.0*x-1.0, PolyN, dPolyN, nDegree)

PolyN  = (-1.0)**(REAL(nDegree))*PolyN
dPolyN = 2.0*(-1.0)**(REAL(nDegree))*dPolyN
        
END SUBROUTINE LegendrePolyNormalized01
!===============================================================================!


!===============================================================================!
SUBROUTINE JacobiRoots(roots)
USE MOD_FiniteVolume_vars,ONLY: nMoms
IMPLICIT NONE
REAL, INTENT(OUT) :: roots(1:nMoms)
INTEGER  :: ii, info
REAL     :: JacobiMatrix(1:nMoms,1:nMoms), bdiag(1:nMoms-1)
INTEGER  :: lwork
REAL, ALLOCATABLE :: work(:) !USE REAL(8)


! Compute recurrence coefficients b_k
DO ii = 2, nMoms
     bdiag(ii-1) = DBLE((ii - 1)*(ii + 1)) / DBLE((2*ii - 1)*(2*ii + 1))
END DO

JacobiMatrix = 0.0
DO ii = 1, nMoms - 1
     JacobiMatrix(ii,   ii+1)  =  sqrt(bdiag(ii))
     JacobiMatrix(ii+1,   ii)  =  sqrt(bdiag(ii))
END DO

! Workspace query for DSYEV
lwork = -1
ALLOCATE(work(1))
CALL DSYEV('N', 'U', nMoms, JacobiMatrix, nMoms, roots, work, lwork, info)
lwork = int(work(1))
DEALLOCATE(work)
ALLOCATE(work(lwork))
CALL DSYEV('N', 'U', nMoms, JacobiMatrix, nMoms, roots, work, lwork, info)
DEALLOCATE(work)

END SUBROUTINE JacobiRoots
!===============================================================================!


!===============================================================================!
END MODULE MOD_PhysicsFrame
!-------------------------------------------------------------------------------!
