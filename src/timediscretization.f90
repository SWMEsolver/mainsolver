!===============================================================================!
MODULE MOD_TimeDiscretization
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE

ABSTRACT INTERFACE
    SUBROUTINE SelectTimeApproximation(t)    
    REAL,INTENT(IN)  :: t
    END SUBROUTINE SelectTimeApproximation
END INTERFACE    

ABSTRACT INTERFACE
    SUBROUTINE SelectTimeStep(deltat)    
    REAL,INTENT(OUT)  :: deltat
    END SUBROUTINE SelectTimeStep
END INTERFACE    

PROCEDURE(SelectTimeApproximation),   pointer :: TimeApproximation => null()
PROCEDURE(SelectTimeStep),            pointer :: TimeStep          => null()

!-------------------------------------------------------------------------------!
INTERFACE TimeDiscretizationByForwardEuler
   MODULE PROCEDURE TimeDiscretizationByForwardEuler
END INTERFACE

INTERFACE TimeDiscretizationBySSPRK4
   MODULE PROCEDURE TimeDiscretizationBySSPRK4
END INTERFACE

INTERFACE TimeDiscretizationByRK65
   MODULE PROCEDURE TimeDiscretizationByRK65
END INTERFACE

INTERFACE TimeDiscretizationByDeC5
   MODULE PROCEDURE TimeDiscretizationByDeC5
END INTERFACE

INTERFACE TimeStepSWE
    MODULE PROCEDURE TimeStepSWE
END INTERFACE

INTERFACE TimeStepGeneralSWME  
  MODULE PROCEDURE TimeStepGeneralSWME 
END INTERFACE

INTERFACE TimeStepHSWME  
  MODULE PROCEDURE TimeStepHSWME 
END INTERFACE

INTERFACE TimeStepSWLME
  MODULE PROCEDURE TimeStepSWLME
END INTERFACE

!-------------------------------------------------------------------------------!
PUBLIC :: TimeDiscretizationByForwardEuler
PUBLIC :: TimeDiscretizationBySSPRK4
PUBLIC :: TimeDiscretizationByRK65
PUBLIC :: TimeDiscretizationByDeC5
PUBLIC :: TimeStepSWE 
PUBLIC :: TimeStepGeneralSWME
PUBLIC :: TimeStepHSWME
PUBLIC :: TimeStepSWLME
! Pointers
PUBLIC :: TimeApproximation
PUBLIC :: TimeStep
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!


!===============================================================================!
SUBROUTINE TimeDiscretizationBySSPRK4(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume,     ONLY: FVSpaceIteration
USE MOD_FiniteVolume_vars,ONLY: U, Ut
USE MOD_FiniteVolume_vars,ONLY: K0, K1, K2, K3, K4, K5
USE MOD_FiniteVolume_vars,ONLY: nXs, nYs, dt
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
REAL            :: tStage
INTEGER         :: ii, jj

! This method can be found in https://doi.org/10.1137/S0036142901389025
! Butcher tableau in Page 489
! Also see: https://arxiv.org/pdf/1708.02595, Page 10

K0 = U(:,1:nXs,1:nYs)
!-------------------------------------------------------------------------------!
CALL FVSpaceIteration(t + 0.0*dt)
K1 = 1.0*K0 + 0.39175222700392*Ut*dt                                   
U(:,1:nXs,1:nYs)  = K1
!-------------------------------------------------------------------------------!
CALL FVSpaceIteration(t + 0.39175222700392*dt)
K2 = 0.44437049406734*K0 + 0.55562950593266*K1 + 0.36841059262959*Ut*dt
U(:,1:nXs,1:nYs) = K2
!-------------------------------------------------------------------------------!
CALL FVSpaceIteration(t + 0.58607968896780*dt)
K3 = 0.62010185138540*K0 + 0.37989814861460*K2 + 0.25189177424738*Ut*dt
U(:,1:nXs,1:nYs) = K3
!-------------------------------------------------------------------------------!
CALL FVSpaceIteration(t + 0.474542364687*dt)
K4 = 0.17807995410773*K0 + 0.82192004589227*K3 + 0.54497475021237*Ut*dt
U(:,1:nXs,1:nYs) = K4
K5 = Ut
!-------------------------------------------------------------------------------!
CALL FVSpaceIteration(t + 0.93501063100924*dt)
U(:,1:nXs,1:nYs) = 0.00683325884039*K0 +  0.51723167208978*K2 + 0.12759831133288*K3 + &
                   0.34833675773694*K4 + (0.08460416338212*K5 + 0.22600748319395*Ut)*dt
!-------------------------------------------------------------------------------!
END SUBROUTINE TimeDiscretizationBySSPRK4
!===============================================================================!


!===============================================================================!
SUBROUTINE TimeDiscretizationByForwardEuler(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume,     ONLY: FVSpaceIteration
USE MOD_FiniteVolume_vars,ONLY: U, Ut, K0
USE MOD_FiniteVolume_vars,ONLY: nXs, nYs, dt
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
REAL            :: tStage
INTEGER         :: ii, jj

K0 = U(:,1:nXs,1:nYs)
CALL FVSpaceIteration(t)
U(:,1:nXs,1:nYs) = K0 + Ut*dt

END SUBROUTINE TimeDiscretizationByForwardEuler
!===============================================================================!


!===============================================================================!
SUBROUTINE TimeDiscretizationByRK65(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume,     ONLY: FVSpaceIteration
USE MOD_FiniteVolume_vars,ONLY: U, Ut, UN0
USE MOD_FiniteVolume_vars,ONLY: K0, K1, K2, K3, K4, K5
USE MOD_FiniteVolume_vars,ONLY: nXs, nYs, dt
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
REAL, DIMENSION(6,5) :: ARK
REAL, DIMENSION(6)   :: bRK, cRK
REAL                 :: tStage
INTEGER              :: ii, jj, nStages

nStages = 6

ARK = reshape((/ 0., 0.25, 0.125, 0.,  0.1875, -3.0/7.0, &
                 0., 0.,   0.125, 0., -0.375,   8.0/7.0, &
                 0., 0.,   0.,   0.5,  0.375,   6.0/7.0, &
                 0., 0.,   0.,   0. , 0.5625, -12.0/7.0, &
                 0., 0.,   0.,   0. , 0.    ,   8.0/7.0 /), shape(ARK))


bRK = (1.0/90.0)*(/ 7.0, 0.0, 32.0, 12.0, 32.0, 7.0/)
cRK = (/ 0.0, 0.25, 0.25, 0.5, 0.75, 1.0 /)
!-------------------------------------------------------------------------------!
UN0 = U(:,1:nXs,1:nYs)
CALL FVSpaceIteration(t + cRK(1)*dt)
K0 = Ut
!-------------------------------------------------------------------------------!
U(:,1:nXs,1:nYs) = UN0 + ARK(2,1)*K0*dt
CALL FVSpaceIteration(t + cRK(2)*dt)
K1 = Ut
!-------------------------------------------------------------------------------!
U(:,1:nXs,1:nYs) = UN0 + (ARK(3,1)*K0 + ARK(3,2)*K1)*dt
CALL FVSpaceIteration(t + cRK(3)*dt)
K2 = Ut
!-------------------------------------------------------------------------------!
U(:,1:nXs,1:nYs) = UN0 + (ARK(4,1)*K0 + ARK(4,2)*K1 + ARK(4,3)*K2)*dt
CALL FVSpaceIteration(t + cRK(4)*dt)
K3 = Ut
!-------------------------------------------------------------------------------!
U(:,1:nXs,1:nYs) = UN0 + (ARK(5,1)*K0 + ARK(5,2)*K1 + ARK(5,3)*K2 + ARK(5,4)*K3)*dt
CALL FVSpaceIteration(t + cRK(5)*dt)
K4 = Ut
!-------------------------------------------------------------------------------!
U(:,1:nXs,1:nYs) = UN0 + (ARK(6,1)*K0 + ARK(6,2)*K1 +  ARK(6,3)*K2 + ARK(6,4)*K3 + ARK(6,5)*K4 )*dt
CALL FVSpaceIteration(t + cRK(6)*dt)
K5 = Ut
!-------------------------------------------------------------------------------!
U(:,1:nXs,1:nYs) = UN0 + (bRK(1)*K0 + bRK(2)*K1 + bRK(3)*K2 + bRK(4)*K3 + bRK(5)*K4 + bRK(6)*K5)*dt

END SUBROUTINE TimeDiscretizationByRK65
!===============================================================================!


!===============================================================================!
SUBROUTINE TimeDiscretizationByDeC5(t) !
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume,     ONLY: FVSpaceIteration
USE MOD_FiniteVolume_vars,ONLY: U, Ut
USE MOD_FiniteVolume_vars,ONLY: FUp, Ua, Up
USE MOD_FiniteVolume_vars,ONLY: nXs, nYs, dt
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
REAL, DIMENSION(4,4) :: thetaDeC
REAL, DIMENSION(4)   :: betaDeC
REAL, DIMENSION(4)   :: tStage
INTEGER, PARAMETER   :: KCorr = 5, MSteps = 4
INTEGER              :: ii, jj, kk

thetaDeC = reshape((/ 0.0,      &
                      0.1103005664791649049, &
                      0.0730327668541684294, &
                      1.0/12.0, &
                      0.0,      &
                      0.1896994335208351257, &
                      0.4505740308958107176, &
                      5.0/12.0, &
                      0.0,      &
                     -0.0339073642291439076, &
                      0.2269672331458314485, &
                      5.0/12.0, &
                      0.0,      &
                      0.0103005664791649201, &
                     -0.0269672331458315692, &
                      1.0/12.0 /), shape(thetaDeC))

! U^m(k) = U^0 + \Delta t \sum theta_{r,m} F(U^r(k-1))
betaDeC = (/ 0.0, 0.2763932022500210639 , 0.7236067977499789361 , 1.0/)
tStage  = t + betaDeC*dt

! To fix sizes because of the ghost cells
DO ii = 1, MSteps
    Ua(ii, :, :, :) = U
    Up(ii, :, :, :) = U
END DO

CALL FVSpaceIteration(tStage(1))

DO ii = 1, MSteps
    FUp(ii,:,:,:) = Ut
END DO
DO kk = 1, KCorr
    IF (kk == KCorr) THEN
        ii = MSteps
        Ua(ii,:,:,:) = Ua(1,:,:,:)
        DO jj = 1, MSteps
            Ua(ii,:,:,:) = Ua(ii,:,:,:) + dt*thetaDeC(ii, jj)*FUp(jj,:,:,:)
        END DO
    ELSE
        DO ii = 2, MSteps
            Ua(ii,:,:,:) = Ua(1,:,:,:)
            DO jj = 1, MSteps
                Ua(ii,:,:,:)  = Ua(ii,:,:,:) + dt*thetaDeC(ii, jj)*FUp(jj,:,:,:)
            END DO
        END DO
    END IF
    IF (kk .NE. KCorr) THEN
        DO ii = 2, MSteps
            U = Ua(ii,:,:,:)
            CALL FVSpaceIteration(tStage(ii))
            FUp(ii,:,:,:) = Ut
        END DO
    END IF
    Up = Ua
END DO

U = Ua(MSteps,:,:,:)

END SUBROUTINE TimeDiscretizationByDeC5
!===============================================================================!


!===============================================================================!
SUBROUTINE TimeStepSWE(dtSWE)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: U
USE MOD_FiniteVolume_vars,ONLY: CFL, MESH_DX, MIN_TIMESTEP, Gravity
USE MOD_FiniteVolume_vars,ONLY: nVar, nDims, nXs, nYs
USE MOD_PhysicsFrame,     ONLY: ConsToPrim
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(OUT) :: dtSWE
!-------------------------------------------------------------------------------!
REAL    :: FastestWaveX, FastestWaveY
REAL    :: LambdaMax(1:2)
REAL    :: Prim(1:nVar), h, vx, vy
INTEGER :: ii, jj, kk

LambdaMax = 0.0
dtSWE     = HUGE(1.0)

#ifdef TWODIM
DO jj = 1,nYs
    DO ii = 1,nXs
        CALL ConsToPrim(U(:,ii,jj), Prim)
        
        h  = Prim(1)
        vx = Prim(2)
        vy = Prim(3)

        FastestWaveX = vx + SQRT(Gravity*h)
        FastestWaveY = vy + SQRT(Gravity*h)

        LambdaMax(1) = MAX(LambdaMax(1),  ABS(FastestWaveX))
        LambdaMax(2) = MAX(LambdaMax(2),  ABS(FastestWaveY))
        
        DO kk = 1,nDims
            dtSWE   = MIN(dtSWE,  MESH_DX(kk)/LambdaMax(kk))
        END DO
    END DO
END DO

dtSWE = CFL*dtSWE

IF (dtSWE .LT. MIN_TIMESTEP) THEN
    dtSWE = MIN_TIMESTEP
END IF

#endif

END SUBROUTINE TimeStepSWE
!===============================================================================!


!===============================================================================!
SUBROUTINE TimeStepGeneralSWME(dtSWME)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: U, CFL, MESH_DX, nVar, nMoms, nXs, nYs
USE MOD_FiniteVolume_vars,ONLY: MIN_TIMESTEP, Gravity
USE MOD_PhysicsFrame,     ONLY: ConsToPrim 
USE MOD_MomentModels
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL, INTENT(OUT) :: dtSWME
!-------------------------------------------------------------------------------!
REAL    :: FastestWaveX, FastestWaveY, wPlus, wMins, maxEV
REAL    :: Prim(1:nVar), sqrWave(1:nXs,nYs)
REAL    :: JA(1:nVar, 1:nVar), JB(1:nVar, 1:nVar)
REAL    :: WR(1:nVar), WI(1:nVar)                ! WR/WI: Real/Imaginary parts of eigenvalues
REAL    :: VL(1:nVar,1:nVar), VR(1:nVar,1:nVar)  ! VL/VR: Left/Right eigenvectors
REAL, ALLOCATABLE :: work(:)
      
INTEGER :: ii, jj, kk, lwork, info
! ABOUT DGEEV: IS FOR DOUBLE PRECISION REAL NUMBERS AS SET IN MAKE
! [OUT]    WR, WI, VL (for left ev), VR, work, info
! [IN,OUT] Matrix A
JA = 0.0; JB = 0.0
WR = 0.0; WI = 0.0
VL = 0.0; VR = 0.0
maxEV = 0.0

dtSWME = HUGE(1.0)
DO ii = 1, nXs
    DO jj = 1, nYs        
        CALL ConsToPrim(U(1:nVar,ii,jj), Prim(1:nVar))                
        !--------------------------------------------------------------------------------
        CALL SystemMatrix(Prim, JA)    
        CALL RemainMatrix(Prim, JB)    
        !--------------------------------------------------------------------------------
        JA = JA + JB
        !--------------------------------------------------------------------------------
        ! STRICT DOUBLE SYNTAX FOR THE DGEEV SOLVER; MODIFICATIONS MAY LEAD TO PROBLEMS
        lwork = -1
        ALLOCATE(work(1))
        CALL DGEEV('N', 'V', nVar, JA, nVar, WR, WI,  VL, nVar, VR, nVar, work, lwork, info)        
        lwork = int(work(1))
        DEALLOCATE(work)
        ALLOCATE(work(lwork))
        CALL DGEEV('N', 'V', nVar, JA, nVar, WR, WI,  VL, nVar, VR, nVar, work, lwork, info)        
        DEALLOCATE(work)
        !--------------------------------------------------------------------------------
        DO kk = 1, nVar
            maxEV = MAX(maxEV, WR(kk))
        END DO                
        dtSWME = MIN(dtSWME, MESH_DX(1)/maxEV)
    END DO
END DO
dtSWME = CFL*dtSWME

END SUBROUTINE TimeStepGeneralSWME
!===============================================================================!


!!===============================================================================!
SUBROUTINE TimeStepHSWME(dtHSWME)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: U, V, CFL, MESH_DX, nVar, nMoms, nXs, nYs
USE MOD_FiniteVolume_vars,ONLY: MIN_TIMESTEP, Gravity, RootsAH
USE MOD_MomentModels
USE MOD_PhysicsFrame,     ONLY: ConsToPrim 
IMPLICIT NONE  
!-------------------------------------------------------------------------------!
REAL, INTENT(OUT) :: dtHSWME
!-------------------------------------------------------------------------------!
REAL    :: wPlus, wMins, maxEV, Prim(1:nVar)
REAL    :: sqrWave(1:nXs,1:nYs), EigValj(1:nMoms+2)      
INTEGER :: ii, jj, kk, lwork, info

!! Paper by Koellermeier & Rominger; doi: 10.4208/cicp.OA-2019-0065; page 17.
maxEV = 0.0;
dtHSWME = HUGE(1.0)

!print*, "heree", SHAPE(V),'shapeeeee', shape(U)
DO ii = 1,nXs
    DO jj = 1,nYs
        CALL ConsToPrim(U(1:nVar,ii,jj), Prim(1:nVar))
                        
        EigValj(1) = Prim(2) + sqrt(Gravity*Prim(1) + Prim(3)*Prim(3))
        EigValj(2) = Prim(2) - sqrt(Gravity*Prim(1) + Prim(3)*Prim(3))

        maxEV = max(maxEV, abs(EigValj(1)))
        maxEV = max(maxEV, abs(EigValj(2)))
        
        DO kk = 1, nMoms           
            EigValj(kk + 2) = Prim(2) + RootsAH(kk)*Prim(3)         
            maxEV = max(maxEV, abs(EigValj(kk + 2))) 
        END DO                
    END DO    
END DO

dtHSWME = CFL*MESH_DX(1)/maxEV

END SUBROUTINE TimeStepHSWME
!===============================================================================!



!!===============================================================================!
SUBROUTINE TimeStepSWLME(dtSWLME)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: U, V, CFL, MESH_DX, nVar, nMoms, nXs, nYs
USE MOD_FiniteVolume_vars,ONLY: MIN_TIMESTEP, Gravity, RootsAH
USE MOD_MomentModels
USE MOD_PhysicsFrame,     ONLY: ConsToPrim 
IMPLICIT NONE  
!-------------------------------------------------------------------------------!
REAL, INTENT(OUT) :: dtSWLME
!-------------------------------------------------------------------------------!
REAL    :: wPlus, wMins, maxEV, Prim(1:nVar)
REAL    :: sqrWave(1:nXs,1:nYs), EigValj(1:nMoms+2)      
INTEGER :: ii, jj, kk, lwork, info

!! Paper by Koellermeier & Rominger; doi: 10.4208/cicp.OA-2019-0065; page 17.
maxEV = 0.0;
dtSWLME = HUGE(1.0)

DO ii = 1,nXs
    DO jj = 1,nYs
        CALL ConsToPrim(U(1:nVar,ii,jj), Prim(1:nVar))
        
        EigValj(1) = Prim(2) + sqrt(Gravity*Prim(1) + Prim(3)*Prim(3))
        EigValj(2) = Prim(2) - sqrt(Gravity*Prim(1) + Prim(3)*Prim(3))

        maxEV = max(maxEV, abs(EigValj(1)))
        maxEV = max(maxEV, abs(EigValj(2)))
        
        DO kk = 1, nMoms           
            EigValj(kk + 2) = Prim(2)
            maxEV = max(maxEV, abs(EigValj(kk + 2))) 
        END DO                
    END DO    
END DO

dtSWLME = CFL*MESH_DX(1)/maxEV

END SUBROUTINE TimeStepSWLME
!===============================================================================!

!===============================================================================!
END MODULE MOD_TimeDiscretization
!-------------------------------------------------------------------------------!
