!===============================================================================!
MODULE MOD_Reconstruction
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ReconstructionX
  MODULE PROCEDURE ReconstructionX
END INTERFACE

INTERFACE ReconstructionY
  MODULE PROCEDURE ReconstructionY
END INTERFACE

!INTERFACE ReconstructionFixX
!  MODULE PROCEDURE ReconstructionFixX
!END INTERFACE

!INTERFACE ReconstructionFixY
!  MODULE PROCEDURE ReconstructionFixY
!END INTERFACE

INTERFACE PositivityLimiterX
  MODULE PROCEDURE PositivityLimiterX
END INTERFACE 

INTERFACE PositivityLimiterY
  MODULE PROCEDURE PositivityLimiterY
END INTERFACE 

INTERFACE MUSCL 
  MODULE PROCEDURE MUSCL 
END INTERFACE

INTERFACE WENO3_SecondSweep 
  MODULE PROCEDURE WENO3_SecondSweep
END INTERFACE

INTERFACE WENO5_SecondSweep 
  MODULE PROCEDURE WENO5_SecondSweep
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: ReconstructionX
PUBLIC :: ReconstructionY
!PUBLIC :: ReconstructionFixX
!PUBLIC :: ReconstructionFixY
PUBLIC :: PositivityLimiterX
PUBLIC :: PositivityLimiterY
PUBLIC :: MUSCL
PUBLIC :: WENO3_SecondSweep 
PUBLIC :: WENO5_SecondSweep 
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!


!===============================================================================!
SUBROUTINE ReconstructionX()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: V, U, WM, WP
USE MOD_FiniteVolume_vars,ONLY: nXs, nYs, nGPsX
USE MOD_FiniteVolume_vars,ONLY: nGhostsX, nGhostsY, MESH_DX
USE MOD_FiniteVolume_vars,ONLY: Ind, ReconstrFlag
IMPLICIT NONE
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage


#ifndef TWODIM
SELECT CASE (ReconstrFlag)
CASE (1)
    DO ii = 0, nXs + 1
        WM(:,1,ii,1) = V(:,ii,1)
        WP(:,1,ii,1) = V(:,ii,1)
    END DO
    
    
CASE (2)

    DO ii = 0, nXs + 1
        CALL MUSCL(V(:,ii-1,1),V(:,ii,1),V(:,ii+1,1),WM(:,1,ii,1),WP(:,1,ii,1),MESH_DX(1))            
    END DO
CASE (3, 4)
    DO ii = 0, nXs + 1
        IF (.NOT. Ind(1, ii, 1)) THEN
            CALL WENO_XDIR( V(:, -nGhostsX + ii:ii + nGhostsX, 1), WM(:,1, ii, 1), WP(:,1, ii, 1), ReconstrFlag)
        END IF
    END DO
CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE (*, *) ErrorMessage
    STOP
END SELECT

#else

SELECT CASE (ReconstrFlag)
CASE (1)
    DO jj = 1, nYs
        DO ii = 0, nXs + 1
            WM(:, nGPsX, ii, jj) = V(:, ii, jj)
            WP(:, nGPsX, ii, jj) = V(:, ii, jj)
        END DO
    END DO    
CASE (2)
    DO jj = 1, nYs
        DO ii = 0, nXs + 1            
            CALL MUSCL(V(:,ii-1,jj), V(:,ii,jj), V(:,ii+1,jj), WM(:,:, ii, jj), WP(:,:, ii, jj), MESH_DX(1))            
        END DO
    END DO
    
CASE (3, 4)
    DO jj = 1, nYs
        DO ii = 0, nXs + 1
            IF (.NOT. Ind(1, ii, jj)) THEN
                CALL WENO_XDIR( V(:, -nGhostsX + ii:ii + nGhostsX, -nGhostsY + jj:jj + nGhostsY), &
                                WM(:,:, ii, jj), WP(:,:, ii, jj), ReconstrFlag)
            END IF
        END DO
    END DO        
CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE (*, *) ErrorMessage
    STOP
END SELECT

#endif


! This subroutine gives the values of WP and WM 
! used in the computation of the flux by the Riemann - Solver
END SUBROUTINE ReconstructionX
!===============================================================================!


!===============================================================================!
SUBROUTINE ReconstructionY()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: V, WM, WP
USE MOD_FiniteVolume_vars,ONLY: nDims, nXs, nYs, nGPsX
USE MOD_FiniteVolume_vars,ONLY: nGhostsX, nGhostsY, MESH_DX
USE MOD_FiniteVolume_vars,ONLY: Ind, ReconstrFlag
IMPLICIT NONE
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: ErrorMessage


SELECT CASE (ReconstrFlag)
CASE (1)
    DO jj = 0, nYs + 1
        DO ii = 1, nXs
            WM(:, nGPsX, ii, jj) = V(:, ii, jj) ! TODO loop over NGPsX?
            WP(:, nGPsX, ii, jj) = V(:, ii, jj)
        END DO
    END DO
    
    
CASE (2)

    DO jj = 0, nYs + 1
        DO ii = 1, nXs
            CALL MUSCL( V(:,ii,jj-1), V(:,ii,jj), V(:,ii,jj+1), WM(:,1, ii, jj), WP(:,1, ii, jj), MESH_DX(nDims))                        
        END DO
    END DO
    
    
CASE (3, 4)
#ifdef TWODIM            
    DO jj = 0, nYs + 1
        DO ii = 1, nXs
            IF (.NOT. Ind(2, ii, jj)) THEN
                CALL WENO_YDIR( V(:, -nGhostsX + ii:ii + nGhostsX, -nGhostsY + jj:jj + nGhostsY), &
                                WM(:,:, ii, jj), WP(:,:, ii, jj), ReconstrFlag)            
            END IF
        END DO
    END DO
#else
    DO ii = 1, nXs
        IF (.NOT. Ind(2, ii, 1)) THEN
            CALL WENO_YDIR( V(:, -nGhostsX + ii:ii + nGhostsX, 1), WM(:,:, ii, 1), WP(:,:, ii, 1), ReconstrFlag)           
        END IF
    END DO
#endif                        
CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE (*, *) ErrorMessage
    STOP
END SELECT

END SUBROUTINE ReconstructionY
!===============================================================================!


!!===============================================================================!
!SUBROUTINE ReconstructionFixX()
!IMPLICIT NONE
!END SUBROUTINE ReconstructionFixX
!!===============================================================================!


!!===============================================================================!
!SUBROUTINE ReconstructionFixY()
!IMPLICIT NONE
!END SUBROUTINE ReconstructionFixY
!!===============================================================================!


!===============================================================================!
SUBROUTINE PositivityLimiterX()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: V, WM, WP
USE MOD_FiniteVolume_vars,ONLY: nGPsX, nXs, nYs
USE MOD_FiniteVolume_vars,ONLY: WeightsGPBnd, wLobatto, Hmin
IMPLICIT NONE   
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP
REAL    :: alpha, beta, csiX, theta, mmin 

DO jj = 1, nYs
    DO ii = 0, nXs
        alpha = 0.
        beta  = 0.
        DO iGP = 1, nGPsX
            alpha = alpha + WeightsGPBnd(iGP)*WM(1, iGP, ii, jj)
            beta  = beta  + WeightsGPBnd(iGP)*WP(1, iGP, ii, jj)
        END DO        
        csiX = (V(1, ii, jj) - wLobatto*alpha - wLobatto*beta)/(1. - 2.*wLobatto)
        DO iGP = 1, nGPsX
            mmin = MIN(csiX, WM(1, iGP, ii, jj), WP(1, iGP, ii, jj))
            IF (V(1, ii, jj) .EQ. mmin) THEN
                theta = 1.
            ELSE
                theta = MIN(1., ABS((V(1, ii, jj) - Hmin)/(V(1, ii, jj) - mmin)))
            END IF
            WP(1, iGP, ii, jj) = V(1, ii, jj) + theta*(WP(1, iGP, ii, jj) - V(1, ii, jj))
            WM(1, iGP, ii, jj) = V(1, ii, jj) + theta*(WM(1, iGP, ii, jj) - V(1, ii, jj))
        END DO

    END DO
END DO

END SUBROUTINE PositivityLimiterX
!===============================================================================!


!===============================================================================!
SUBROUTINE PositivityLimiterY()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: V, WM, WP
USE MOD_FiniteVolume_vars,ONLY: nGPsY, nXs, nYs
USE MOD_FiniteVolume_vars,ONLY: WeightsGPBnd, wLobatto, Hmin
IMPLICIT NONE 
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, jGP
REAL    :: alpha, beta, csiY, theta, mmin 

DO jj = 0, nYs
    DO ii = 1, nXs
        alpha = 0.0; beta = 0.0        
        DO jGP = 1, nGPsY
            alpha = alpha + WeightsGPBnd(jGP)*WM(1, jGP, ii, jj)
            beta  = beta  + WeightsGPBnd(jGP)*WP(1, jGP, ii, jj)
        END DO
        csiY = (V(1, ii, jj) - wLobatto*alpha - wLobatto*beta)/(1.-2.*wLobatto)
        DO jGP = 1, nGPsY
            mmin = MIN(csiY, WM(1, jGP, ii, jj), WP(1, jGP, ii, jj))

            IF (V(1, ii, jj) .EQ. mmin) THEN
                theta = 1.
            ELSE
                theta = MIN(1., ABS((V(1, ii, jj) - Hmin)/(V(1, ii, jj) - mmin)))
            END IF
            WP(1, jGP, ii, jj) = V(1, ii, jj) + theta*(WP(1, jGP, ii, jj) - V(1, ii, jj))
            WM(1, jGP, ii, jj) = V(1, ii, jj) + theta*(WM(1, jGP, ii, jj) - V(1, ii, jj))
        END DO

    END DO
END DO

END SUBROUTINE PositivityLimiterY
!===============================================================================!


!===============================================================================!
SUBROUTINE MUSCL(PrimM, Prim0, PrimP, WWM, WWP,dx)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: dx
REAL,INTENT(IN)  :: PrimM(1:nVar), Prim0(1:nVar), PrimP(1:nVar)
REAL,INTENT(OUT) :: WWM(1:nVar)
REAL,INTENT(OUT) :: WWP(1:nVar)
!-------------------------------------------------------------------------------!
REAL             :: SM, SP, SLOPE, DuM, DuP, Ratio
INTEGER          :: iVar

SM = 0.0; SP = 0.0
! Here nGPsX = 1 and nGhostsX = 1
DuM = 0.0; DuP = 0.0
Ratio = 0.0 + 0.0*dx

DO iVar = 1, nVar
    DuM = Prim0(iVar) - PrimM(iVar) 
    DuP = PrimP(iVar) - Prim0(iVar) 
    
    ! Factor out DuP:
    !     = MAX(0.0, MIN(1.0, Ratio))  -->  DuM/DuP = sgn(DuP)*DuM/ABS(DuP)
    SLOPE = MAX(0.0, MIN(ABS(DuP), SIGN(1.0,DuP)*DuM))*SIGN(1.0,DuP)
    ! Other Slope
    !SLOPE = MAX(0.0, MIN(2.0*DuP, MIN(2.0*DuM,   0.5*(DuP + DuM))))
   
    ! Note that DuP has been cancelled out with the factor from DuP
    WWM(iVar) = Prim0(iVar) + 0.5*SLOPE
    WWP(iVar) = Prim0(iVar) - 0.5*SLOPE
END DO
        
END SUBROUTINE MUSCL
!===============================================================================!


!===============================================================================!
SUBROUTINE WENO_XDIR(Vst,WWM,WWP,WhichReconstruction)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nGPsX, nGhostsX, nGhostsY, WENOEXP, WENOEPS
IMPLICIT NONE
!-------------------------------------------------------------------------------!
#ifdef TWODIM
    REAL,INTENT(IN)    :: Vst(1:nVar,-nGhostsX:nGhostsX,-nGhostsY:nGhostsY)
    REAL,INTENT(OUT)   :: WWM(1:nVar,1:nGPsX)
    REAL,INTENT(OUT)   :: WWP(1:nVar,1:nGPsX)
#else
    REAL,INTENT(IN)    :: Vst(1:nVar,-nGhostsX:nGhostsX)
    REAL,INTENT(OUT)   :: WWM(1:nVar)
    REAL,INTENT(OUT)   :: WWP(1:nVar)
#endif
INTEGER,INTENT(IN) :: WhichReconstruction
!-------------------------------------------------------------------------------!
REAL               :: beta0, beta1, gamma0, gamma1, denominator, omega0, omega1
REAL               :: VtempM(1:nVar,-nGhostsX:nGhostsX)
REAL               :: VtempP(1:nVar,-nGhostsX:nGhostsX)
INTEGER            :: iVar, jj
CHARACTER(LEN=255) :: ErrorMessage


! Vst is the solution restricted at the stencil 
! so that (i - Ghosts: i + Ghosts) becomes simply -Ghosts:Ghosts

! Better to reconstruct the conservative variables
! Check in the PVM scheme

SELECT CASE (WhichReconstruction)
CASE (3)
    DO iVar = 1, nVar        
#ifdef TWODIM
        DO jj = -nGhostsY, nGhostsY
            CALL WENO3_FirstSweep(Vst(iVar,:, jj), VtempM(iVar, jj), VtempP(iVar, jj))
        END DO             
        CALL WENO3_SecondSweep(VtempM(iVar,:), WWM(iVar, :)) ! Minus, values at the quadrature point along the edge in y
        CALL WENO3_SecondSweep(VtempP(iVar,:), WWP(iVar, :)) ! Plus,  values at the quadrature point along the edge in y
#else        
        beta0  = (Vst(iVar,0)  - Vst(iVar,1))**2
        beta1  = (Vst(iVar,-1) - Vst(iVar,0))**2
        gamma0 = 1.0/3.0
        gamma1 = 2.0/3.0
        denominator = gamma0*(WENOEPS + beta0)**WENOEXP + gamma1*(WENOEPS + beta1)**WENOEXP
        omega0 = gamma0*(WENOEPS +  beta0)**(WENOEXP)/denominator
        omega1 = gamma1*(WENOEPS +  beta1)**(WENOEXP)/denominator
        
        !----------------------------------------------------------------
        ! cell edge x_{j-1/2}  value
        !----------------------------
        ! u_{j-1/2}^+
        WWM(iVar) = .5*(omega0 * (1.0*Vst(iVar,0) + Vst(iVar,-1))) + & 
                    .5*(omega1 * (3.0*Vst(iVar,0) - Vst(iVar,1)))

        !----------------------------------------------------------------
        ! cell edge x_{j+1/2}  value          
        !---------------------------- 
        ! u_{j+1/2}^-       
        WWP(iVar) = .5*(omega0 * (3.0*Vst(iVar,0) - Vst(iVar,-1))) + & 
                    .5*(omega1 * (1.0*Vst(iVar,0) + Vst(iVar,1)))
        ! Be aware that here I have swapped WWP with WWM from the order given originally in the WENO3_first       
        ! One has to investigate why this is happening, there might be an issue with indexes or with a sign somewhere.
#endif
    END DO

CASE (4)
    DO iVar = 1, nVar
#ifdef TWODIM
        DO jj = -nGhostsY, nGhostsY
            CALL WENO5_FirstSweep( Vst(iVar,:, jj), VtempM(iVar, jj), VtempP(iVar, jj))
        END DO
        CALL WENO5_SecondSweep(VtempM(iVar,:), WWM(iVar, :))
        CALL WENO5_SecondSweep(VtempP(iVar,:), WWP(iVar, :))
#else  
        CALL WENO5_FirstSweep( Vst(iVar,:), WWM(iVar), WWP(iVar))    
#endif
    END DO
CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE (*, *) ErrorMessage
    STOP
END SELECT

END SUBROUTINE WENO_XDIR
!===============================================================================!


!===============================================================================!
SUBROUTINE WENO_YDIR(V,WWM,WWP,WhichReconstruction)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nGPsX, nGhostsX, nGhostsY
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: V(1:nVar,-nGhostsX:nGhostsX,-nGhostsY:nGhostsY)
REAL,INTENT(OUT)   :: WWM(1:nVar,1:nGPsX)
REAL,INTENT(OUT)   :: WWP(1:nVar,1:nGPsX)
INTEGER,INTENT(IN) :: WhichReconstruction
!-------------------------------------------------------------------------------!
REAL               :: VtempM(1:nVar,-nGhostsX:nGhostsX)
REAL               :: VtempP(1:nVar,-nGhostsX:nGhostsX)
INTEGER            :: iVar, ii
CHARACTER(LEN=255) :: ErrorMessage

SELECT CASE (WhichReconstruction)
CASE (3)
    DO iVar = 1, nVar
        DO ii = -nGhostsX, nGhostsX
            CALL WENO3_FirstSweep(V(iVar, ii, -nGhostsY:nGhostsY), VtempM(iVar, ii), VtempP(iVar, ii))
        END DO
        CALL WENO3_SecondSweep(VtempM(iVar, -nGhostsX:nGhostsX), WWM(iVar, :))
        CALL WENO3_SecondSweep(VtempP(iVar, -nGhostsX:nGhostsX), WWP(iVar, :))
    END DO
CASE (4)
    DO iVar = 1, nVar
        DO ii = -nGhostsX, nGhostsX
            CALL WENO5_FirstSweep( V(iVar, ii, -nGhostsY:nGhostsY), VtempM(iVar, ii), VtempP(iVar, ii))
        END DO
        CALL WENO5_SecondSweep(VtempM(iVar, -nGhostsX:nGhostsX), WWM(iVar, :))
        CALL WENO5_SecondSweep(VtempP(iVar, -nGhostsX:nGhostsX), WWP(iVar, :))
    END DO
CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE (*, *) ErrorMessage
    STOP
END SELECT

END SUBROUTINE WENO_YDIR
!===============================================================================!


!===============================================================================!
SUBROUTINE WENO3_FirstSweep(Q,WWM,WWP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nGhostsX, WENOEPS, WENOEXP
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhostsX:nGhostsX)
REAL,INTENT(OUT) :: WWM
REAL,INTENT(OUT) :: WWP
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2
REAL             :: beta1, beta2
REAL             :: gamma1, gamma2
REAL             :: omega1, omega2
REAL             :: W1, W2

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (Q(-1) - Q(0))**2.0
beta2 = (Q(0)  - Q(1))**2.0

!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!
! Linear Weights
gamma1 = 2.0/3.0
gamma2 = 1.0/3.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

W1 = 0.5*(   Q(-1) + Q(0))
W2 = 0.5*(3.0*Q(0) - Q(1))
WWM = omega1*W1 + omega2*W2

!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!

! Linear Weights
gamma1 = 1.0/3.0
gamma2 = 2.0/3.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

! Reconstructed Polynomial
W1 = 0.5*(-Q(-1) + 3.0*Q(0))
W2 = 0.5*( Q(0) +     Q(1))
WWP = omega1*W1 + omega2*W2

END SUBROUTINE WENO3_FirstSweep
!===============================================================================!


!===============================================================================!
SUBROUTINE WENO3_SecondSweep(Q,W)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nGhostsX, nGPsX, WENOEPS, WENOEXP
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhostsX:nGhostsX)
REAL,INTENT(OUT) :: W(1:nGPsX)
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2
REAL             :: beta1, beta2
REAL             :: gamma1, gamma2
REAL             :: omega1, omega2
REAL             :: W1, W2

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (Q(-1) - Q(0))**2.0
beta2 = (Q(0)  - Q(1))**2.0

!------------------------------!
! Point: x_{j-1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
gamma1 = 0.5
gamma2 = 0.5

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

! Reconstructed Polynomial
W1   = (1.0/6.0)*(SQRT(3.0)*Q(-1) + 6.0*Q(0) - SQRT(3.0)*Q(0))
W2   = (1.0/6.0)*(SQRT(3.0)*Q(0)  + 6.0*Q(0) - SQRT(3.0)*Q(1))
W(1) = omega1*W1 + omega2*W2

!------------------------------!
! Point: x_{j+1/(2*sqrt(3))}   !
!------------------------------!

! Linear Weights
gamma1 = 0.5
gamma2 = 0.5

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2)
omega2 = alpha2/(alpha1 + alpha2)

! Reconstructed Polynomial
W1   = (1.0/6.0)*(-SQRT(3.0)*Q(-1) + 6.0*Q(0) + SQRT(3.0)*Q(0))
W2   = (1.0/6.0)*(-SQRT(3.0)*Q(0)  + 6.0*Q(0) + SQRT(3.0)*Q(1))
W(2) = omega1*W1 + omega2*W2

END SUBROUTINE WENO3_SecondSweep
!===============================================================================!


!===============================================================================!
SUBROUTINE WENO5_FirstSweep(Q,WWM,WWP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nGhostsX, WENOEPS, WENOEXP
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhostsX:nGhostsX)
REAL,INTENT(OUT) :: WWM
REAL,INTENT(OUT) :: WWP
!-------------------------------------------------------------------------------!
REAL             :: alpha1, alpha2, alpha3
REAL             :: beta1,  beta2,  beta3
REAL             :: gamma1, gamma2, gamma3
REAL             :: omega1, omega2, omega3
REAL             :: W1, W2, W3

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
                 + 11.0*Q(-2)*Q(0) - 31.0*Q(-1)*Q(0) + 10.0*Q(0)*Q(0))
beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(0) + 13.0*Q(0)*Q(0) &
                 +  5.0*Q(-1)*Q(1) - 13.0*Q(0)*Q(1) +  4.0*Q(1)*Q(1))
beta3 = (1.0/3.0)*(10.0*Q(0)*Q(0) - 31.0*Q(0)*Q(1) + 25.0*Q(1)*Q(1) &
                 + 11.0*Q(0)*Q(+2) - 19.0*Q(1)*Q(+2) +  4.0*Q(+2)*Q(+2))

!------------------------------!
! WM: x_{i-1/2}                !
!------------------------------!

! Linear Weights
gamma1 = 3.0/10.0
gamma2 = 3.0/5.0
gamma3 = 1.0/10.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

W1 = (1.0/6.0)*(    -Q(-2) + 5.0*Q(-1) + 2.0*Q(0))
W2 = (1.0/6.0)*( 2.0*Q(-1) + 5.0*Q(0) -     Q(1))
W3 = (1.0/6.0)*(11.0*Q(0) - 7.0*Q(1) + 2.0*Q(+2))
WWM = omega1*W1 + omega2*W2 + omega3*W3

!------------------------------!
! WP: x_{i+1/2}                !
!------------------------------!

! Linear Weights
gamma1 = 1.0/10.0
gamma2 = 3.0/5.0
gamma3 = 3.0/10.0

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

W1 = (1.0/6.0)*(2.0*Q(-2) - 7.0*Q(-1) + 11.0*Q(0))
W2 = (1.0/6.0)*(   -Q(-1) + 5.0*Q(0) +  2.0*Q(1))
W3 = (1.0/6.0)*(2.0*Q(0)  + 5.0*Q(1) -      Q(+2))
WWP = omega1*W1 + omega2*W2 + omega3*W3

END SUBROUTINE WENO5_FirstSweep
!===============================================================================!


!===============================================================================!
SUBROUTINE WENO5_SecondSweep(Q,W) !4nGPS
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nGhostsX, nGPsX, WENOEPS, WENOEXP
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Q(-nGhostsX:nGhostsX)
REAL,INTENT(OUT) :: W(1:nGPsX)
!-------------------------------------------------------------------------------!
REAL :: alpha1, alpha2, alpha3
REAL :: beta1,  beta2,  beta3
REAL :: gamma1, gamma2, gamma3
REAL :: omega1, omega2, omega3
REAL :: W1, W2, W3

!------------------------------!
! Common Smoothness Indicators !
!------------------------------!
beta1 = (1.0/3.0)*( 4.0*Q(-2)*Q(-2) - 19.0*Q(-2)*Q(-1) + 25.0*Q(-1)*Q(-1) &
                 + 11.0*Q(-2)*Q(0) - 31.0*Q(-1)*Q(0) + 10.0*Q(0)*Q(0))
beta2 = (1.0/3.0)*( 4.0*Q(-1)*Q(-1) - 13.0*Q(-1)*Q(0) + 13.0*Q(0)*Q(0) &
                 +  5.0*Q(-1)*Q(1) - 13.0*Q(0)*Q(1) +  4.0*Q(1)*Q(1))
beta3 = (1.0/3.0)*(10.0*Q(0)*Q(0) - 31.0*Q(0)*Q(1) + 25.0*Q(1)*Q(1) &
                 + 11.0*Q(0)*Q(+2) - 19.0*Q(1)*Q(+2) +  4.0*Q(+2)*Q(+2))

!--------------------------------------------!
! Point: x_{j-1/2*sqrt(3/7+2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 = 0.2658420974778319
gamma2 = 0.6112504900322486
gamma3 = 0.1229074124899195

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1 = -0.1642562761719537*Q(-2)+0.7590807081409336*Q(-1)+0.4051755680310201*Q(0)+0.0000000000000000*Q(1)+0.0000000000000000*Q(2)
W2 = +0.0000000000000000*Q(-2)+0.2663118796250726*Q(-1)+0.8979443965468811*Q(0)-0.1642562761719537*Q(1)+0.0000000000000000*Q(2)
W3 = +0.0000000000000000*Q(-2)+0.0000000000000000*Q(-1)+1.6968800354220990*Q(0)-0.9631919150471715*Q(1)+0.2663118796250726*Q(2)
W(1) = (omega1*W1 + omega2*W2 + omega3*W3)

!--------------------------------------------!
! Point: x_{j-1/2*sqrt(3/7-2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights

gamma1 = 0.1281641584355059
gamma2 = 0.5219691498503563
gamma3 = 0.3498666917141379

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)


W1 = -0.1122135388132497*Q(-2)+0.3944175994189276*Q(-1)+0.7177959393943222*Q(0)+0.0000000000000000*Q(1)+0.0000000000000000*Q(2)
W2 = +0.0000000000000000*Q(-2)+0.0577769829791784*Q(-1)+1.0544365558340714*Q(0)-0.1122135388132497*Q(1)+0.0000000000000000*Q(2)
W3 = +0.0000000000000000*Q(-2)+0.0000000000000000*Q(-1)+1.2277675047716066*Q(0)-0.2855444877507849*Q(1)+0.0577769829791784*Q(2)

W(2) = (omega1*W1 + omega2*W2 + omega3*W3)


!--------------------------------------------!
! Point: x_{j+1/2*sqrt(3/7-2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 = 0.3498666917141379
gamma2 = 0.5219691498503563
gamma3 = 0.1281641584355059

alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial

W1 = +0.0577769829791784*Q(-2)-0.2855444877507849*Q(-1)+1.2277675047716066*Q(0)+0.0000000000000000*Q(1)+0.0000000000000000*Q(2)
W2 = +0.0000000000000000*Q(-2)-0.1122135388132497*Q(-1)+1.0544365558340714*Q(0)+0.0577769829791784*Q(1)+0.0000000000000000*Q(2)
W3 = +0.0000000000000000*Q(-2)+0.0000000000000000*Q(-1)+0.7177959393943222*Q(0)+0.3944175994189276*Q(1)-0.1122135388132497*Q(2)

W(3) = (omega1*W1 + omega2*W2 + omega3*W3)

!--------------------------------------------!
! Point: x_{j+1/2*sqrt(3/7+2/7*sqrt(6/5))}   !
!--------------------------------------------!

! Linear Weights
gamma1 = 0.1229074124899195
gamma2 = 0.6112504900322486
gamma3 = 0.2658420974778319


alpha1 = gamma1/(WENOEPS + beta1)**WENOEXP
alpha2 = gamma2/(WENOEPS + beta2)**WENOEXP
alpha3 = gamma3/(WENOEPS + beta3)**WENOEXP

! Nonlinear Weights
omega1 = alpha1/(alpha1 + alpha2 + alpha3)
omega2 = alpha2/(alpha1 + alpha2 + alpha3)
omega3 = alpha3/(alpha1 + alpha2 + alpha3)

! Reconstructed Polynomial
W1 = +0.2663118796250726*Q(-2)-0.9631919150471715*Q(-1)+1.6968800354220990*Q(0)+0.0000000000000000*Q(1)+0.0000000000000000*Q(2)
W2 = +0.0000000000000000*Q(-2)-0.1642562761719537*Q(-1)+0.8979443965468811*Q(0)+0.2663118796250726*Q(1)+0.0000000000000000*Q(2)
W3 = +0.0000000000000000*Q(-2)+0.0000000000000000*Q(-1)+0.4051755680310201*Q(0)+0.7590807081409336*Q(1)-0.1642562761719537*Q(2)
W(4) = (omega1*W1 + omega2*W2 + omega3*W3)

END SUBROUTINE WENO5_SecondSweep
!===============================================================================!


!===============================================================================!
END MODULE MOD_Reconstruction
!-------------------------------------------------------------------------------!
