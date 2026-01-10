!===============================================================================!
MODULE MOD_FiniteVolume
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE

ABSTRACT INTERFACE
    SUBROUTINE SelectViscousMatrix(MatA, MatB)
    USE MOD_FiniteVolume_vars,ONLY: nVar
    REAL,INTENT(IN)  :: MatA(1:nVar, 1:nVar)
    REAL,INTENT(OUT) :: MatB(1:nVar, 1:nVar)        
    END SUBROUTINE SelectViscousMatrix
END INTERFACE    

ABSTRACT INTERFACE
    SUBROUTINE SelectIntegralPath(Sval, PrimL, PrimR, Phi, dPhi)
    USE MOD_FiniteVolume_vars,ONLY: nVar
    REAL,INTENT(IN)  :: Sval, PrimL(1:nVar), PrimR(1:nVar)
    REAL,INTENT(OUT) :: Phi(1:nVar), dPhi
    END SUBROUTINE SelectIntegralPath
END INTERFACE    


ABSTRACT INTERFACE
    SUBROUTINE SelectNonConsMethod(PrimL,PrimR, DvalM, DvalP)
    USE MOD_FiniteVolume_vars, ONLY: nVar
    REAL,INTENT(IN)  :: PrimL(1:nVar)
    REAL,INTENT(IN)  :: PrimR(1:nVar)
    REAL,INTENT(OUT) :: DvalM(1:nVar)
    REAL,INTENT(OUT) :: DvalP(1:nVar)
    END SUBROUTINE SelectNonConsMethod
END INTERFACE

ABSTRACT INTERFACE
    subroutine SelectConservativeFlux(PrimL, PrimR, NormVect, TangVect, Flux)
    USE MOD_FiniteVolume_vars,ONLY: nVar, nGPsX    
    REAL,INTENT(IN)  :: PrimL(1:nVar,1:nGPsX)
    REAL,INTENT(IN)  :: PrimR(1:nVar,1:nGPsX)
    REAL,INTENT(IN)  :: NormVect(1:2,1:nGPsX)
    REAL,INTENT(IN)  :: TangVect(1:2,1:nGPsX)  
    REAL,INTENT(OUT) :: Flux(1:nVar,1:nGPsX)
    END SUBROUTINE SelectConservativeFlux
END INTERFACE 
PROCEDURE(SelectViscousMatrix),   pointer :: ViscousMatrix    => null()
PROCEDURE(SelectIntegralPath),    pointer :: PathIntegral     => null()
PROCEDURE(SelectConservativeFlux),pointer :: ConservativeFlux => null()
PROCEDURE(SelectNonConsMethod),   pointer :: NonConsMethod    => null()
!-------------------------------------------------------------------------------!
INTERFACE FVSpaceIteration
    MODULE PROCEDURE FVSpaceIteration
END INTERFACE

INTERFACE FillInitialConditions
    MODULE PROCEDURE FillInitialConditions
END INTERFACE

INTERFACE InitializeWBVariables
    MODULE PROCEDURE InitializeWBVariables
END INTERFACE

INTERFACE InitializeFiniteVolume
    MODULE PROCEDURE InitializeFiniteVolume
END INTERFACE

INTERFACE FinalizeFiniteVolume
    MODULE PROCEDURE FinalizeFiniteVolume
END INTERFACE

INTERFACE ZeroSolver  
  MODULE PROCEDURE ZeroSolver 
END INTERFACE

INTERFACE ViscousMatrixLxF 
  MODULE PROCEDURE ViscousMatrixLxF 
END INTERFACE

INTERFACE ViscousMatrixLxW  
  MODULE PROCEDURE ViscousMatrixLxW 
END INTERFACE

INTERFACE ViscousMatrixForce  
  MODULE PROCEDURE ViscousMatrixForce 
END INTERFACE

INTERFACE ViscousMatrixHLL
  MODULE PROCEDURE ViscousMatrixHLL
END INTERFACE

INTERFACE LinearPath
  MODULE PROCEDURE LinearPath
END INTERFACE

INTERFACE QuadraticPath
  MODULE PROCEDURE QuadraticPath
END INTERFACE

INTERFACE PowerlawMinusPath
  MODULE PROCEDURE PowerlawMinusPath
END INTERFACE

INTERFACE PowerlawPlusPath
  MODULE PROCEDURE PowerlawPlusPath
END INTERFACE

INTERFACE RiemannSolver
    MODULE PROCEDURE RiemannSolver
END INTERFACE

INTERFACE PVM
    MODULE PROCEDURE PVM
END INTERFACE

INTERFACE EigenValueMethod
    MODULE PROCEDURE EigenValueMethod
END INTERFACE

INTERFACE ZeroNonCons 
  MODULE PROCEDURE ZeroNonCons 
END INTERFACE

!-------------------------------------------------------------------------------!
PUBLIC :: FVSpaceIteration
PUBLIC :: FillInitialConditions
PUBLIC :: InitializeWBVariables
PUBLIC :: InitializeFiniteVolume
PUBLIC :: FinalizeFiniteVolume
! Specific subroutines
PUBLIC :: ZeroSolver
PUBLIC :: ViscousMatrixLxF
PUBLIC :: ViscousMatrixLxW
PUBLIC :: ViscousMatrixForce
PUBLIC :: ViscousMatrixHLL
PUBLIC :: LinearPath
PUBLIC :: QuadraticPath
PUBLIC :: PowerlawMinusPath
PUBLIC :: PowerlawPlusPath
PUBLIC :: RiemannSolver
PUBLIC :: PVM
PUBLIC :: EigenValueMethod
PUBLIC :: ZeroNonCons
! Pointers:
PUBLIC :: ViscousMatrix
PUBLIC :: PathIntegral
PUBLIC :: ConservativeFlux
PUBLIC :: NonConsMethod
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!


!===============================================================================!
SUBROUTINE FVSpaceIteration(t)
!-------------------------------------------------------------------------------!
USE MOD_PhysicsFrame,   ONLY: BoundaryConditions, PrimToCons 
USE MOD_Reconstruction, ONLY: ReconstructionX, ReconstructionY
USE MOD_Reconstruction, ONLY: PositivityLimiterX, PositivityLimiterY
USE MOD_ShocksIndicator,ONLY: ShocksIndicatorX, ShocksIndicatorY

USE MOD_FiniteVolume_vars,ONLY: U, V, Ut, SSFriction, SWB, WM, WP, FX, FY
USE MOD_FiniteVolume_vars,ONLY: AuxFlux, FXWB, FYWB
USE MOD_FiniteVolume_vars,ONLY: nVar, nXs, nYs, MESH_DX, nGPsX, nGPsY, IndexPass
USE MOD_FiniteVolume_vars,ONLY: NormVectX, TangVectX, NormVectY, TangVectY
USE MOD_FiniteVolume_vars,ONLY: WeightsGPBnd
USE MOD_FiniteVolume_vars,ONLY: ModelFlag, NCFluxFlag, SplittingFlag

USE MOD_FiniteVolume_vars,ONLY: AuxDplus, AuxDmins
USE MOD_FiniteVolume_vars,ONLY: DplusX, DminsX, DplusY, DminsY, Bath, BvecX, BvecXrot

USE MOD_PhysicsFrame,     ONLY: ConsToPrim

USE MOD_MomentModels
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: t
!-------------------------------------------------------------------------------!
REAL    :: NNv(1:2,1:nGPsX)
REAL    :: TTv(1:2,1:nGPsX)
REAL    :: BathDminsX(1:nVar), AuxBathDminsX(1:nVar,1:nGPsX,1:nXs,1:nYs)
REAL    :: BathDplusX(1:nVar), AuxBathDplusX(1:nVar,1:nGPsX,1:nXs,1:nYs), Dummy(1:nVar)
REAL    :: Vleft(1:nVar), Vright(1:nVar)
INTEGER :: ii, jj, iGP, jGP

SSFriction = 0.0
FX = 0.0; FY = 0.0; AuxFlux = 0.0;
AuxDplus = 0.0; DplusX = 0.0;  DminsX = 0.0
AuxDmins = 0.0; DplusY = 0.0;  DminsY = 0.0

BathDminsX = 0.0; AuxBathDminsX = 0.0
BathDplusX = 0.0; AuxBathDplusX = 0.0

Vleft  = 0.0
Vright = 0.0
Dummy  = 0.0
Ut     = 0.0

!!-------------------------------------------------------------------------------!
!! Update V (needs to be here when having multiple steps time integration):
DO jj = 1, nYs
    DO ii = 1, nXs
        CALL ConsToPrim( U(:, ii, jj),  V(:, ii, jj) )
    END DO
END DO
!!-------------------------------------------------------------------------------!
CALL BoundaryConditions(t)
!CALL ShocksIndicatorX()    
CALL ReconstructionX()    
!CALL PositivityLimiterX()  

! Both the higher order weno and muscl schemes need to be fixed/adapted
! For now run with reconstruction flag 1
#ifndef TWODIM
DO ii = 0, nXs        
    NNv = NormVectX(:,:,ii,1)
    TTv = TangVectX(:,:,ii,1)
    VLeft  = WP(:,1,ii,1);    
    VRight = WM(:,1,ii+1,1);
    CALL NonConsMethod( VLeft,  VRight, DminsX(:,ii,1), DplusX(:,ii,1))           
END DO
DO ii = 1, nXs                    
    ! Non-conservative flux: DP(u_{i-1},u_i) + DM(u_i,u_{i+1})        
    Ut(:,ii,1) = Ut(:,ii,1) - (FX(:,ii,1) - FX(:,ii-1,1))/Mesh_DX(1)
    Ut(:,ii,1) = Ut(:,ii,1) - (DplusX(:,ii-1,1) + DminsX(:,ii,1))/Mesh_DX(1)
END DO
!------------------------------------------------------------------------------------------
! Two dimensional case
#else

DO jj = 1, nYs
    DO ii = 0, nXs        
        NNv = NormVectX(:,:,ii,jj)
        TTv = TangVectX(:,:,ii,jj)

        CALL ConservativeFlux(  WP(:,:,ii,jj),  WM(:,:,ii+1,jj),  NNv,  TTv,  AuxFlux)
        
        DO iGP = 1, nGPsX
            
            VLeft  = WP(:,iGP,ii,jj);
            VRight = WM(:,iGP+1,ii,jj);
            CALL NonConsMethod( VLeft,  VRight, AuxDmins, AuxDplus)
            FX(:,ii,jj) = FX(:,ii,jj) + weightsGPBnd(iGP)*AuxFlux(:, iGP)            
            
            DplusX(:,ii,jj) = DplusX(:,ii,jj) + weightsGPBnd(iGP)*AuxDplus
            DminsX(:,ii,jj) = DminsX(:,ii,jj) + weightsGPBnd(iGP)*AuxDmins      
              
        END DO
    END DO
END DO

CALL ShocksIndicatorY()
CALL ReconstructionY()
CALL PositivityLimiterY()

DO jj = 0, nYs
    DO ii = 1, nXs
    
        NNv = NormVectY(:,:,ii,jj)
        TTv = TangVectY(:,:,ii,jj)        
        CALL ConservativeFlux(WP(:,:,ii,jj), WM(:,:,ii,jj+1), NNv, TTv, AuxFlux)
        
        DO jGP = 1, nGPsY
            CALL NonConsMethod(WP(:,jGP,ii,jj), WM(:,jGP,ii,jj+1), AuxDmins, AuxDplus)                                  
            FY(:,ii,jj)     =  FY(:,ii,jj)      +  weightsGPBnd(jGP)*AuxFlux(:, jGP)            
            DplusY(:,ii,jj) =  DplusY(:,ii,jj)  +  weightsGPBnd(jGP)*AuxDplus
            DminsY(:,ii,jj) =  DminsY(:,ii,jj)  +  weightsGPBnd(jGP)*AuxDmins
        END DO
        
        
    END DO
END DO

!-------------------------------------------------------------------------------!
!CALL SourceFuncSWE2D(t) ! Gives S ! To include for Bathymetry
!#ifdef WELLBALANCED
!    FX = FX - FXWB
!    FY = FY - FYWB
!    S  = S  - SWB
!#endif

DO jj = 1, nYs
    DO ii = 1, nXs                    
                 
        ! X direction
        Ut(:,ii,jj) = Ut(:,ii,jj) - (FX(:,ii,jj) - FX(:,ii-1,jj))/Mesh_DX(1)
        Ut(:,ii,jj) = Ut(:,ii,jj) - (DplusX(:,ii-1,jj) + DminsX(:,ii,jj))/Mesh_DX(1)
        ! Y direction
        Ut(:,ii,jj) = Ut(:,ii,jj) - (FY(:,ii,jj) - FY(:,ii,jj-1))/Mesh_DX(2)
        Ut(:,ii,jj) = Ut(:,ii,jj) - (DplusY(:,ii,jj-1) + DminsY(:,ii,jj))/Mesh_DX(2)
        
    END DO
END DO
#endif

END SUBROUTINE FVSpaceIteration
!===============================================================================!


!===============================================================================!
SUBROUTINE PVM(PrimL,PrimR, DvalM, DvalP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars, ONLY: nVar
USE MOD_PhysicsFrame,      ONLY: PrimToCons
USE MOD_FiniteVolume_vars, ONLY: QuadPoints, QuadWeights, nGPath
USE MOD_MomentModels 
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: PrimL(1:nVar)
REAL,INTENT(IN)  :: PrimR(1:nVar)
REAL,INTENT(OUT) :: DvalM(1:nVar)
REAL,INTENT(OUT) :: DvalP(1:nVar)
!-------------------------------------------------------------------------------!
REAL             :: ConsL(1:nVar), ConsR(1:nVar)
REAL             :: matrixA(1:nVar, 1:nVar),  matrixAtemp(1:nVar, 1:nVar)
REAL             :: matrixB(1:nVar, 1:nVar)
REAL             :: matrixQ(1:nVar, 1:nVar)
REAL             :: Phi(1:nVar), dPhi
INTEGER          :: iGPk
!------------------------------------------------------------------------------------------------!

ConsL   = 0.0
ConsR   = 0.0
matrixA = 0.0
matrixB = 0.0
matrixQ = 0.0
matrixAtemp = 0.0

DvalM = 0.0
DvalP = 0.0
!-------------------------------------------------------------------------------!  
 
CALL PrimToCons(PrimL, ConsL)
CALL PrimToCons(PrimR, ConsR)  
!-------------------------------------------------------------------------------
! Non conservative fluctuation; Linear path; middle point quadrature:

DO iGPk = 1, nGPath

    CALL PathIntegral(QuadPoints(iGPk), PrimL, PrimR, Phi, dPhi)

    CALL SystemMatrix(Phi, matrixAtemp)
    CALL RemainMatrix(Phi, matrixB    )

    !-------------------------------------------------------------------------------
    matrixAtemp = matrixAtemp + matrixB
    !-------------------------------------------------------------------------------
    matrixA = matrixA + QuadWeights(iGPk) * matrixAtemp*dPhi
END DO

! matrixA here is what is called the linearization of the Roe matrix
CALL ViscousMatrix(matrixA, matrixQ)

DvalM = 0.5*MATMUL( matrixA - matrixQ, ConsR - ConsL) 
DvalP = 0.5*MATMUL( matrixA + matrixQ, ConsR - ConsL)

    
END SUBROUTINE PVM
!===============================================================================!


!===============================================================================!
SUBROUTINE ViscousMatrixLxF(matrixA, matrixQ)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, dt, Identity, MESH_DX
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: matrixA(1:nVar, 1:nVar)
REAL,INTENT(OUT) :: matrixQ(1:nVar, 1:nVar)

matrixQ = 0.0
matrixQ = (Mesh_DX(1)/dt)*Identity    

END SUBROUTINE ViscousMatrixLxF
!===============================================================================!


!===============================================================================!
SUBROUTINE ViscousMatrixLxW(matrixA, matrixQ)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, dt, MESH_DX
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: matrixA(1:nVar, 1:nVar)
REAL,INTENT(OUT) :: matrixQ(1:nVar, 1:nVar)

matrixQ = 0.0
matrixQ = (dt/Mesh_DX(1))*MATMUL(matrixA, matrixA)
   
END SUBROUTINE ViscousMatrixLxW
!===============================================================================!


!===============================================================================!
SUBROUTINE ViscousMatrixForce(matrixA, matrixQ)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, dt, Identity, MESH_DX
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: matrixA(1:nVar, 1:nVar)
REAL,INTENT(OUT) :: matrixQ(1:nVar, 1:nVar)

matrixQ = 0.0
matrixQ = .5*(Mesh_DX(1)/dt)*Identity + .5*(dt/Mesh_DX(1))*MATMUL(matrixA, matrixA)       

END SUBROUTINE ViscousMatrixForce
!===============================================================================!


!===============================================================================!
SUBROUTINE ViscousMatrixHLL(matrixA, matrixQ)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, dt, Identity, MESH_DX
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: matrixA(1:nVar, 1:nVar)
REAL,INTENT(OUT) :: matrixQ(1:nVar, 1:nVar)
!-------------------------------------------------------------------------------!
REAL    :: maxEV, minEV, aa0, aaN, JA(1:nVar,1:nVar)
REAL    :: VL(1:nVar,1:nVar), VR(1:nVar,1:nVar) ! VL/VR: Left/Right eigenvectors
REAL    :: WR(1:nVar), WI(1:nVar)               ! WR/WI: Real/Imaginary parts of eigenvalues     
INTEGER :: kk, lwork, info
REAL, ALLOCATABLE :: work(:)

matrixQ = 0.0
maxEV = 0.0
minEV = 0.0

WR = 0.0; WI = 0.0
VL = 0.0; VR = 0.0

JA = matrixA
!-------------------------------------------------------------------------------!
lwork = -1
ALLOCATE(work(1))
CALL DGEEV('N', 'V', nVar, JA, nVar, WR, WI, VL, nVar, VR, nVar, work, lwork, info)        
lwork = int(work(1))
DEALLOCATE(work)
ALLOCATE(work(lwork))
CALL DGEEV('N', 'V', nVar, JA, nVar, WR, WI, VL, nVar, VR, nVar, work, lwork, info)        
DEALLOCATE(work)
!-------------------------------------------------------------------------------!
maxEV = WR(1)
minEV = WR(1)
DO kk = 1, nVar  
    maxEV = MAX(maxEV, WR(kk))
    minEV = MIN(minEV, WR(kk))
END DO        
aa0 = (maxEV*ABS(minEV) - minEV*ABS(maxEV))/(maxEV - minEV)
aaN = (ABS(maxEV) - ABS(minEV))/(maxEV - minEV)

matrixQ = aa0*Identity + aaN*matrixA       


END SUBROUTINE ViscousMatrixHLL
!===============================================================================!


!===============================================================================!
SUBROUTINE EigenValueMethod(PrimL,PrimR, DvalM, DvalP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars, ONLY: nVar
USE MOD_PhysicsFrame,      ONLY: PrimToCons
USE MOD_FiniteVolume_vars, ONLY: QuadPoints, QuadWeights, nGPath
USE MOD_MomentModels 
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: PrimL(1:nVar)
REAL,INTENT(IN)  :: PrimR(1:nVar)
REAL,INTENT(OUT) :: DvalM(1:nVar)
REAL,INTENT(OUT) :: DvalP(1:nVar)
!-------------------------------------------------------------------------------!
REAL             :: matrixAplus(1:nVar, 1:nVar)
REAL             :: matrixAmins(1:nVar, 1:nVar)
REAL             :: ConsL(1:nVar), ConsR(1:nVar)
REAL             :: matrixA(1:nVar, 1:nVar),  matrixAtemp(1:nVar, 1:nVar)
REAL             :: matrixB(1:nVar, 1:nVar)
REAL             :: matrixQ(1:nVar, 1:nVar)
REAL             :: Phi(1:nVar), dPhi
REAL    :: DiagEVplus(1:nVar,1:nVar), DiagEVmins(1:nVar,1:nVar), maxEV, minEV
REAL    :: VL(1:nVar,1:nVar), VR(1:nVar,1:nVar) ! VL/VR: Left/Right eigenvectors
REAL    :: VRinv(1:nVar,1:nVar) 
REAL    :: WR(1:nVar), WI(1:nVar)               ! WR/WI: Real/Imaginary parts of eigenvalues     
INTEGER :: kk, lwork, info,  iGPk, ii, ipiv(1:nVar)
REAL, ALLOCATABLE :: work(:)
!---------------------------------------------------------------------------------------------------------------------
REAL :: h_m, u_m, Prim_m(1:nVar), D, SL, SR, alpha0, alpha1, auxH, auxV(1:nVar), auxAC(1:nVar), S(1:nVar)
!REAL, ALLOCATABLE :: A_inv(:,:), ipiv(:)

!-------------------------------------------------------------------------------!
WR = 0.0; WI = 0.0
VL = 0.0; VR = 0.0
matrixA = 0.0
matrixB = 0.0
matrixQ = 0.0
matrixAtemp = 0.0
DiagEVplus  = 0.0
DiagEVmins  = 0.0
DvalM = 0.0
DvalP = 0.0
!-------------------------------------------------------------------------------!   
CALL PrimToCons(PrimL, ConsL)
CALL PrimToCons(PrimR, ConsR)
!-------------------------------------------------------------------------------!
DO iGPk = 1, nGPath

    CALL PathIntegral(QuadPoints(iGPk), PrimL, PrimR, Phi, dPhi)

    CALL SystemMatrix(Phi, matrixAtemp)
    CALL RemainMatrix(Phi, matrixB    )

    !---------------------------------------------------------------------------!
    matrixAtemp = matrixAtemp + matrixB
    !---------------------------------------------------------------------------!
    matrixA = matrixA + QuadWeights(iGPk) * matrixAtemp*dPhi
END DO
matrixAtemp = matrixA
!-------------------------------------------------------------------------------!

lwork = -1
ALLOCATE(work(1))
CALL DGEEV('N', 'V', nVar, matrixAtemp, nVar, WR, WI, VL, nVar, VR, nVar, work, lwork, info)        
lwork = int(work(1))
DEALLOCATE(work)
ALLOCATE(work(lwork))
CALL DGEEV('N', 'V', nVar, matrixAtemp, nVar, WR, WI, VL, nVar, VR, nVar, work, lwork, info)        
DEALLOCATE(work)
!-------------------------------------------------------------------------------!
maxEV = WR(1)
minEV = WR(1)
DO kk = 1, nVar
    DiagEVplus(kk,kk) = max(0.0,WR(kk))
    DiagEVmins(kk,kk) = min(0.0,WR(kk))
    maxEV = max(maxEV, WR(kk))
    minEV = max(minEV, WR(kk))
END DO        
VRinv = VR

CALL dgetrf(nVar, nVar, VRinv, nVar, ipiv, info)

! --- Step 2: Query optimal workspace size ---
lwork = -1
allocate(work(1))
call dgetri(nVar, VRinv, nVar, ipiv, work, lwork, info)
lwork = int(work(1))
deallocate(work)
allocate(work(lwork))

! --- Step 3: Compute matrix inverse ---
call dgetri(nVar, VRinv, nVar, ipiv, work, lwork, info)
DEALLOCATE(work)

matrixAplus = MATMUL(VR, MATMUL(DiagEVplus, VRinv))
matrixAmins = MATMUL(VR, MATMUL(DiagEVmins, VRinv))

DvalM = MATMUL( matrixAmins, ConsR - ConsL) 
DvalP = MATMUL( matrixAplus, ConsR - ConsL)
!---------------------------------------------------------------
END SUBROUTINE EigenValueMethod
!===============================================================================!


!===============================================================================!
SUBROUTINE LinearPath(Sval, PrimL, PrimR, Phi, dPhi)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Sval, PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT) :: Phi(1:nVar), dPhi

Phi  = PrimL + Sval*(PrimR - PrimL)
dPhi = 1.0 ! Without the factor R-L

END SUBROUTINE LinearPath
!===============================================================================!


!===============================================================================!
SUBROUTINE QuadraticPath(Sval, PrimL, PrimR, Phi, dPhi)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Sval, PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT) :: Phi(1:nVar), dPhi

Phi  = PrimL + Sval*Sval*(PrimR - PrimL)
dPhi = 2*Sval

END SUBROUTINE QuadraticPath
!===============================================================================!


!===============================================================================!
SUBROUTINE PowerlawMinusPath(Sval, PrimL, PrimR, Phi, dPhi)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Sval, PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT) :: Phi(1:nVar), dPhi
REAL :: PowerParam
PowerParam = 0.7

Phi  = PrimL + (Sval**(PowerParam))*(PrimR - PrimL)
dPhi = PowerParam*Sval**(PowerParam - 1.0)

END SUBROUTINE PowerlawMinusPath
!===============================================================================!


!===============================================================================!
SUBROUTINE PowerlawPlusPath(Sval, PrimL, PrimR, Phi, dPhi)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Sval, PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT) :: Phi(1:nVar), dPhi
REAL :: PowerParam
PowerParam = 0.7

Phi  = PrimR + ((1.0 - Sval)**(PowerParam))*(PrimL - PrimR)
dPhi = PowerParam*(1.0-Sval)**(PowerParam - 1.0)

END SUBROUTINE PowerlawPlusPath
!===============================================================================!


!===============================================================================!
SUBROUTINE BathySource(BathL, BathR, PrimL, PrimR, BathDvalM, BathDvalP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, EpsilonRatio, Gravity, AngleIncl
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: PrimL(1:nVar), BathL
REAL,INTENT(IN)  :: PrimR(1:nVar), BathR
REAL,INTENT(OUT) :: BathDvalM(1:nVar)
REAL,INTENT(OUT) :: BathDvalP(1:nVar)
!-------------------------------------------------------------------------------!
REAL :: AuxVector(1:nVar)

BathDvalM = 0.0 
BathDvalP = 0.0
AuxVector = 0.0; AuxVector(2) = 1.0

!-----------------------------------------------------------------------------
! To write a general version with the quadrature rule and viscosity matrix
BathDvalM = .5* EpsilonRatio* COS(AngleIncl) * (PrimR(1) + PrimL(1))*(BathR - BathL)*AuxVector
BathDvalP = .5* EpsilonRatio* COS(AngleIncl) * (PrimR(1) + PrimL(1))*(BathR - BathL)*AuxVector

END SUBROUTINE BathySource
!===============================================================================!



!===============================================================================!
SUBROUTINE FillInitialConditions()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: U, V
USE MOD_FiniteVolume_vars,ONLY: nVar, nXs, nYs, MESH_X0, MESH_X1, MESH_DX, AngleIncl
USE MOD_FiniteVolume_vars,ONLY: nGPsX, nGPsY, MeshGP, WeightsGP, Bath, dBath, BvecX, BvecXrot
USE MOD_FiniteVolume_vars,ONLY: InitialFlag, BathymetryFlag
USE MOD_PhysicsFrame,ONLY: ConsToPrim, ExactFunction, Bathymetry, DerivativeBathymetry
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL, DIMENSION(nVar,  nGPsX, nGPsY) :: Utemp
REAL, DIMENSION(nGPsX, nGPsY)        :: Btemp
REAL, DIMENSION(nGPsX, nGPsY)        :: dBtemp
REAL :: AuxdBathL, AuxdBathR

INTEGER :: ii, jj, iGP, jGP

!Not affected by the moments
U = 0.  ! Cons variable
V = 0.  ! Prim variable
Utemp = 0.
Btemp = 0.

! Take care of this in 1D
DO jj = 1, nYs
    DO ii = 1, nXs
        ! Loops on quadrature points, x and y axes:
        DO iGP = 1, nGPsX
            DO jGP = 1, nGPsY
                ! compute cell average of conservative variables
                CALL ExactFunction(0.0, MeshGP(:, ii, jj, iGP, jGP), Utemp(:, iGP, jGP))
                U(:, ii, jj) = U(:, ii, jj) + WeightsGP(iGP, jGP)*Utemp(:, iGP, jGP)

                ! compute cell average of bathymetry
                Btemp(iGP, jGP) = Bathymetry(MeshGP(:, ii, jj, iGP, jGP))
                Bath(ii, jj)    = Bath(ii, jj) + WeightsGP(iGP, jGP)*Btemp(iGP, jGP)
                
                ! compute cell average of derivative of bathymetry        
                dBtemp(iGP, jGP) = DerivativeBathymetry(MeshGP(:, ii, jj, iGP, jGP))
                dBath(ii, jj)    = dBath(ii, jj) + WeightsGP(iGP, jGP)*dBtemp(iGP, jGP)       
                
            END DO
        END DO
        
        CALL ConsToPrim(U(:, ii, jj), V(:, ii, jj))
    END DO
    
    
END DO

! Take care of this in 1D
IF (BathymetryFlag .EQ. 4) THEN
    DO ii = 1, nXs        
        AuxdBathL = DerivativeBathymetry( (/Mesh_X0(1) + REAL(ii-1) * MESH_DX(1), 0.0 /))                                   
        AuxdBathR = DerivativeBathymetry( (/Mesh_X0(1) + REAL(ii) * MESH_DX(1),   0.0 /))                                           
        dBath(ii, 1) = .5*(AuxdBathL + AuxdBathR)                                  
    END DO
    CALL ConsToPrim(U(:, ii, jj), V(:, ii, jj))
END IF


!-------------------------------------------------------------------------------!
! For now only valid in 1D
BvecX(0) = Bathymetry((/ MESH_X0(1)/))
DO ii = 1, nXs
    BvecX(ii) = Bathymetry((/MESH_X0(1) + (REAL(ii) - 0.5)*MESH_DX(1) /)) 
END DO
BvecX(nXs+1) = Bathymetry( (/MESH_X1(1)/))
!-------------------------------------------------------------------------------!

BvecXrot(0) = BvecX(0) - MESH_X0(1)*TAN(AngleIncl)
DO ii = 1, nXs
    BvecXrot(ii) = BvecX(ii) - (MESH_X0(1) + (REAL(ii) - 0.5)*MESH_DX(1))*TAN(AngleIncl)
END DO
BvecXrot(nXs+1) = BvecX(nXs+1) - MESH_X1(1)*TAN(AngleIncl)

END SUBROUTINE FillInitialConditions
!===============================================================================!


!===============================================================================!
SUBROUTINE InitializeWBVariables()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: U, V, UtWB, Ut
USE MOD_FiniteVolume_vars,ONLY: FXWB, FX, FYWB, FY, SWB, SSFriction
USE MOD_FiniteVolume_vars,ONLY: nVar, nXs, nYs, nGPsX, nGPsY
USE MOD_FiniteVolume_vars,ONLY: InitialFlag
USE MOD_FiniteVolume_vars,ONLY: MeshGP, WeightsGP
USE MOD_PhysicsFrame,ONLY: ExactFunctionWB
USE MOD_PhysicsFrame,ONLY: ConsToPrim
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL, DIMENSION(nVar, nGPsX, nGPsY) :: Utemp
INTEGER :: ii, jj, iGP, jGP

! compute UWB
Utemp = 0.
DO jj = 1, nYs
    DO ii = 1, nXs
        DO iGP = 1, nGPsX
            DO jGP = 1, nGPsY
                CALL ExactFunctionWB(MeshGP(:, ii, jj, iGP, jGP), Utemp(:, iGP, jGP))
                U(:, ii, jj) = U(:, ii, jj) + WeightsGP(iGP, jGP)*Utemp(:, iGP, jGP)
            END DO
        END DO
        CALL ConsToPrim(U(:, ii, jj), V(:, ii, jj))
    END DO
END DO

CALL FVSpaceIteration(0.)

UtWB = Ut
FXWB = FX
FYWB = FY
SWB  = SSFriction

END SUBROUTINE InitializeWBVariables
!===============================================================================!


!===============================================================================!
SUBROUTINE InitializeFiniteVolume()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nGPsX, nGPsY, nDims, nXs, nYs, nMoms
USE MOD_FiniteVolume_vars,ONLY: MeshNodes, MeshBary, MeshGP, Bath, dBath, BvecX, BvecXrot
USE MOD_FiniteVolume_vars,ONLY: NormVectX, NormVectY, TangVectX, TangVectY
USE MOD_FiniteVolume_vars,ONLY: U, V, Ut, WM, WP, Up, Ua
USE MOD_FiniteVolume_vars,ONLY: FX, FY, AuxFlux, FUp, SSFriction
USE MOD_FiniteVolume_vars,ONLY: UN0, K0, K1, K2, K3, K4, K5
USE MOD_FiniteVolume_vars,ONLY: WeightsGP, WeightsGPBnd, nGhostsX, nGhostsY, Ind
USE MOD_FiniteVolume_vars,ONLY: FXWB, FYWB, UtWB,  SWB 

USE MOD_FiniteVolume_vars,ONLY: PrimRefState1, PrimRefState2, MESH_DX
USE MOD_FiniteVolume_vars,ONLY: PrimRefState3, PrimRefState4

USE MOD_FiniteVolume_vars,ONLY: AuxDplus, AuxDmins
USE MOD_FiniteVolume_vars,ONLY: DplusX, DminsX, DplusY, DminsY
USE MOD_FiniteVolume_vars,ONLY: NNv, TTv, RootsAH
!--------------------------------------------------------------------------------------! 
IMPLICIT NONE
!--------------------------------------------------------------------------------------! 
INTEGER            :: ii, jj, iSparseVect, iRowStart, K, L, PathFlag 
CHARACTER(LEN=255) :: ErrorMessage

ALLOCATE(PrimRefState1(1:nVar),   PrimRefState2(1:nVar),                       & 
         PrimRefState3(1:nVar),   PrimRefState4(1:nVar), MESH_DX(1:nDims))

!--------------------------------------------------------------------------------------!
ALLOCATE(MeshNodes(1:nDims, 0:nXs, 0:nYs), & 
         MeshBary( 1:nDims, 1:nXs, 1:nYs))
ALLOCATE(MeshGP(1:nDims,  1:nXs, 1:nYs, 1:nGPsX, 1:nGPsY))
!--------------------------------------------------------------------------------------!
ALLOCATE(WeightsGP(1:nGPsX, 1:nGPsY),  WeightsGPBnd(1:nGPsX))
!--------------------------------------------------------------------------------------!
ALLOCATE(NormVectX(1:2,  1:nGPsX,  0:nXs,  1:nYs))
ALLOCATE(TangVectX(1:2,  1:nGPsX,  0:nXs,  1:nYs))
ALLOCATE(NormVectY(1:2,  1:nGPsX,  1:nXs,  0:nYs))
ALLOCATE(TangVectY(1:2,  1:nGPsX,  1:nXs,  0:nYs))
!--------------------------------------------------------------------------------------!
ALLOCATE(U(1:nVar, -nGhostsX:nXs+nGhostsX+1,  -nGhostsY:nYs+nGhostsY+1), &
         V(1:nVar, -nGhostsX:nXs+nGhostsX+1,  -nGhostsY:nYs+nGhostsY+1) )
!--------------------------------------------------------------------------------------!
ALLOCATE( Ut(1:nVar,  1:nXs,  1:nYs))
ALLOCATE( Bath(1:nXs,  1:nYs), dBath(1:nXs,  1:nYs), BvecX(0:nXs+1), BvecXrot(0:nXs+1) )
ALLOCATE( WM(1:nVar, 1:nGPsX, 0:nXs+1, 0:nYs+1), & 
          WP(1:nVar, 1:nGPsX, 0:nXs+1, 0:nYs+1) )
!--------------------------------------------------------------------------------------!
! Source:
ALLOCATE(SSFriction(1:nVar, 1:nXs, 1:nYs))
!--------------------------------------------------------------------------------------!
! Fluxes:
ALLOCATE(FX(1:nVar, 0:nXs, 1:nYs),  FY(1:nVar, 1:nXs, 0:nYs))
ALLOCATE(AuxFlux(1:nVar, 1:nGPsX))
!--------------------------------------------------------------------------------------!
ALLOCATE(Ind(1:2, 0:nXs + 1, 0:nYs + 1))
ALLOCATE(UN0(1:nVar, 1:nXs, 1:nYs),   K0(1:nVar, 1:nXs, 1:nYs), &
          K1(1:nVar, 1:nXs, 1:nYs),   K2(1:nVar, 1:nXs, 1:nYs), &
          K3(1:nVar, 1:nXs, 1:nYs),   K4(1:nVar, 1:nXs, 1:nYs), &
          K5(1:nVar, 1:nXs, 1:nYs) )
ALLOCATE(Ua(1:4,1:nVar,1:nXs,1:nYs),   Up(1:4,1:nVar,1:nXs,1:nYs))
ALLOCATE(FUp(1:4, 1:nVar, 1:nXs, 1:nYs))
!--------------------------------------------------------------------------------------!
ALLOCATE(AuxDplus(1:nVar), AuxDmins(1:nVar),  &
         DplusX(1:nVar, 0:nXs, 1:nYs), DminsX(1:nVar, 0:nXs, 1:nYs), &
         DplusY(1:nVar, 1:nXs, 0:nYs), DminsY(1:nVar, 1:nXs, 0:nYs))
         
ALLOCATE(NNv(1:2,1:nGPsX), TTv(1:2,1:nGPsX) )
!--------------------------------------------------------------------------------------!


#ifdef WELLBALANCED
ALLOCATE(UtWB(1:nVar, 1:nXs, 1:nYs),  SWB(1:nVar, 1:nXs, 1:nYs)  &
         FXWB(1:nVar, 0:nXs, 1:nYs), FYWB(1:nVar, 1:nXs, 0:nYs))
#endif
!--------------------------------------------------------------------------------------!
Ind = .FALSE.
WM   = 0.0;  WP = 0.0;  AuxFlux = 0.0;
U    = 0.0;  V  = 0.0;  SSFriction = 0.0;      Ut = 0.0;  FX = 0.0;  FY = 0.0
UN0  = 0.0;  K0 = 0.0;  K1 = 0.0;  K2  = 0.0;  K3 = 0.0;  K4 = 0.0;  K5 = 0.0
Bath = 0.0;  Ua = 0.0;  Up = 0.0;  FUp = 0.0
dBath = 0.0; BvecX = 0.0; BvecXrot = 0.0

AuxDplus = 0.0; AuxDmins = 0.0;
DplusX = 0.0; DminsX = 0.0; DplusY = 0.0; DminsY = 0.0;
!--------------------------------------------------------------------------------------!
#ifdef WELLBALANCED
UtWB = 0.0;  SWB = 0.0;  FXWB = 0.0;  FYWB = 0.0
#endif

END SUBROUTINE InitializeFiniteVolume
!===============================================================================!


!===============================================================================!
SUBROUTINE FinalizeFiniteVolume()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: MeshNodes, MeshBary, MeshGP, Bath, dBath, BvecX, BvecXrot
USE MOD_FiniteVolume_vars,ONLY: WeightsGP, WeightsGPBnd, Ind
USE MOD_FiniteVolume_vars,ONLY: NormVectX, NormVectY, TangVectX, TangVectY
USE MOD_FiniteVolume_vars,ONLY: U, V, Ut, WM, WP, Ua, Up
USE MOD_FiniteVolume_vars,ONLY: FX, FY, AuxFlux, FUp, SSFriction
USE MOD_FiniteVolume_vars,ONLY: UN0, K0, K1, K2, K3, K4, K5

USE MOD_FiniteVolume_vars,ONLY: PrimRefState1, PrimRefState2, MESH_DX
USE MOD_FiniteVolume_vars,ONLY: PrimRefState3, PrimRefState4, VarNameVisu
USE MOD_FiniteVolume_vars,ONLY: Identity, SubIdentity, VecDj, VectorBeta, Ccoef

USE MOD_FiniteVolume_vars,ONLY: AuxDplus, AuxDmins
USE MOD_FiniteVolume_vars,ONLY: DplusX, DminsX, DplusY, DminsY
USE MOD_FiniteVolume_vars,ONLY: NNv, TTv, RootsAH

USE MOD_FiniteVolume_vars,ONLY: QuadSourceWeights, QuadSourceNodes, LegenPoly, DxLegenPoly
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: UtWB, FYWB, FXWB, SWB
!-------------------------------------------------------------------------------!
IMPLICIT NONE
DEALLOCATE(MeshNodes, MeshBary, MeshGP, Bath, dBath, BvecX, BvecXrot)
DEALLOCATE(WeightsGP, WeightsGPBnd, Ind)
DEALLOCATE(NormVectX, NormVectY, TangVectX, TangVectY)
DEALLOCATE(U, V, Ut, WM, WP, Ua, Up)
DEALLOCATE(FX, FY, AuxFlux, FUp, SSFriction)
DEALLOCATE(UN0, K0, K1, K2, K3, K4, K5)

DEALLOCATE(PrimRefState1, PrimRefState2, PrimRefState3, PrimRefState4, &
           MESH_DX, Identity, SubIdentity, VecDj, VectorBeta, VarNameVisu)

DEALLOCATE(Ccoef)
           
DEALLOCATE(AuxDplus, AuxDmins, DplusX, DminsX, DplusY, DminsY)
DEALLOCATE(NNv, TTv)
DEALLOCATE(RootsAH)

DEALLOCATE(QuadSourceWeights, QuadSourceNodes, LegenPoly, DxLegenPoly)

#ifdef WELLBALANCED
DEALLOCATE(UtWB, FXWB, FYWB, SWB)
#endif

END SUBROUTINE FinalizeFiniteVolume
!===============================================================================!


!===============================================================================!
SUBROUTINE ZeroSolver(PrimL, PrimR, NormVect, TangVect, Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nGPsX, nDims
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: PrimL(1:nVar,1:nGPsX)
REAL,INTENT(IN)    :: PrimR(1:nVar,1:nGPsX)
REAL,INTENT(IN)    :: NormVect(1:2,1:nGPsX)
REAL,INTENT(IN)    :: TangVect(1:2,1:nGPsX)
REAL,INTENT(OUT)   :: Flux(1:nVar,1:nGPsX)
Flux = 0.0
END SUBROUTINE ZeroSolver
!===============================================================================!


!===============================================================================!
SUBROUTINE ZeroNonCons(PrimL, PrimR, DvalM, DvalP)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: PrimL(1:nVar)
REAL,INTENT(IN)  :: PrimR(1:nVar)
REAL,INTENT(OUT) :: DvalM(1:nVar)
REAL,INTENT(OUT) :: DvalP(1:nVar)
DvalM = 0.0; DvalP = 0.0
END SUBROUTINE ZeroNonCons
!===============================================================================!


!===============================================================================!
SUBROUTINE RiemannSolver(PrimL, PrimR, NormVect, TangVect, Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars, ONLY: nVar, nGPsX, nDims, Gravity
USE MOD_PhysicsFrame,      ONLY: PrimToCons
USE MOD_MomentModels,      ONLY: FluxSWE1D
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: PrimL(1:nVar,1:nGPsX)
REAL,INTENT(IN)    :: PrimR(1:nVar,1:nGPsX)
REAL,INTENT(IN)    :: NormVect(1:2,1:nGPsX)
REAL,INTENT(IN)    :: TangVect(1:2,1:nGPsX)
REAL,INTENT(OUT)   :: Flux(1:nVar,1:nGPsX)
!-------------------------------------------------------------------------------!
REAL             :: PrimLL(1:nVar,1:nGPsX), PrimRR(1:nVar,1:nGPsX)
REAL             :: ConsLL(1:nVar,1:nGPsX), ConsRR(1:nVar,1:nGPsX)
INTEGER          :: iGP
REAL             :: FluxL(1:nVar), FluxR(1:nVar)
REAL             :: LambdaMax, fastestL, fastestR

! TODO Modify the projection thinking in the moments from 4, nVar by pairs

#ifdef TWODIM
DO iGP = 1, nGPsX
    ! Rotating the vector quantities       !
    ! Change the (u,v) to (n,t)            !
    !---------------------------------------
    ! Left:
    PrimLL(1, iGP) = PrimL(1, iGP)
    PrimLL(2, iGP) = NormVect(1, iGP)*PrimL(2, iGP) + NormVect(2, iGP)*PrimL(3, iGP)
    PrimLL(3, iGP) = TangVect(1, iGP)*PrimL(2, iGP) + TangVect(2, iGP)*PrimL(3, iGP)
    !-------------------------------------------------------------------------------
    ! Right:
    PrimRR(1, iGP) = PrimR(1, iGP)
    PrimRR(2, iGP) = NormVect(1, iGP)*PrimR(2, iGP) + NormVect(2, iGP)*PrimR(3, iGP)
    PrimRR(3, iGP) = TangVect(1, iGP)*PrimR(2, iGP) + TangVect(2, iGP)*PrimR(3, iGP)
    !------------------------------------------------------------------------------

    CALL PrimToCons(PrimLL(:, iGP), ConsLL(:, iGP))
    CALL PrimToCons(PrimRR(:, iGP), ConsRR(:, iGP))

    ! Rusanov's flux:
    CALL FluxSWE1D(PrimLL(:, iGP), FluxL)
    CALL FluxSWE1D(PrimRR(:, iGP), FluxR)
    
    fastestL = PrimLL(2,iGP) + SQRT(Gravity*PrimLL(1,iGP))    
    fastestR = PrimRR(2,iGP) + SQRT(Gravity*PrimRR(1,iGP))

    LambdaMax = MAX(ABS(fastestL), ABS(fastestR))

    Flux(:, iGP) = 0.5*((FluxL + FluxR) - LambdaMax*(ConsRR(:, iGP) - ConsLL(:, iGP)))

    ! Rotating back the vector quantities  !
    Flux(2:3, iGP) = NormVect(:, iGP) * Flux(2, iGP) + TangVect(:, iGP) * Flux(3, iGP)
END DO
#endif

END SUBROUTINE RiemannSolver
!===============================================================================!


!===============================================================================!
END MODULE MOD_FiniteVolume
!-------------------------------------------------------------------------------!
