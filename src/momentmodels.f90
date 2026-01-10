!===============================================================================!
MODULE MOD_MomentModels
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
! Define interface for the procedures to be pointed to
ABSTRACT INTERFACE
    SUBROUTINE SelectSystemMatrix(Prim, MatA)
    USE MOD_FiniteVolume_vars,ONLY: nVar
    REAL,INTENT(IN)  :: Prim(1:nVar)
    REAL,INTENT(OUT) :: MatA(1:nVar, 1:nVar)
    END SUBROUTINE SelectSystemMatrix
    
    SUBROUTINE SelectRemainMatrix(Prim, MatB)
    USE MOD_FiniteVolume_vars,ONLY: nVar
    REAL,INTENT(IN)  :: Prim(1:nVar)
    REAL,INTENT(OUT) :: MatB(1:nVar, 1:nVar)        
    END SUBROUTINE SelectRemainMatrix    

    SUBROUTINE SelectFlux1D(Prim, Flux)
    USE MOD_FiniteVolume_vars,ONLY: nVar
    REAL,INTENT(IN)  :: Prim(1:nVar)
    REAL,INTENT(OUT) :: Flux(1:nVar)    
    END SUBROUTINE SelectFlux1D

    SUBROUTINE SelectSourceTerm(Prim, SourceSW)
    USE MOD_FiniteVolume_vars,ONLY: nVar
    REAL,INTENT(IN)  :: Prim(1:nVar)
    REAL,INTENT(OUT) :: SourceSW(1:nVar)
    END SUBROUTINE SelectSourceTerm

    SUBROUTINE SelectSourceSplit(Prim, SourceSW)
    USE MOD_FiniteVolume_vars,ONLY: nVar
    REAL,INTENT(IN)  :: Prim(1:nVar)
    REAL,INTENT(OUT) :: SourceSW(1:nVar)
    END SUBROUTINE SelectSourceSplit

    !---------------------------------------------------------------------------!    
    SUBROUTINE SelectFlux2DX(Prim, Flux)
    USE MOD_FiniteVolume_vars,ONLY: nVar
    REAL,INTENT(IN)  :: Prim(1:nVar)
    REAL,INTENT(OUT) :: Flux(1:nVar)        
    END SUBROUTINE SelectFlux2DX
        
    SUBROUTINE SelectFlux2DY(Prim, Flux)
    USE MOD_FiniteVolume_vars,ONLY: nVar
    REAL,INTENT(IN)  :: Prim(1:nVar)
    REAL,INTENT(OUT) :: Flux(1:nVar)       
    END SUBROUTINE SelectFlux2DY    
END INTERFACE
!-------------------------------------------------------------------------------!
PROCEDURE(SelectSystemMatrix), pointer :: SystemMatrix => null()
PROCEDURE(SelectRemainMatrix), pointer :: RemainMatrix => null()
PROCEDURE(SelectSourceTerm),   pointer :: SourceTerm   => null()
PROCEDURE(SelectSourceSplit),  pointer :: SourceSplit  => null()
PROCEDURE(SelectFlux1D),       pointer :: Flux1D  => null()
PROCEDURE(SelectFlux2DX),      pointer :: Flux2DX => null()
PROCEDURE(SelectFlux2DY),      pointer :: Flux2DY => null()
!-------------------------------------------------------------------------------!
INTERFACE JacobianSWME1D
    MODULE PROCEDURE JacobianSWME1D
END INTERFACE

INTERFACE JacobianHSWME1D
    MODULE PROCEDURE JacobianHSWME1D
END INTERFACE

INTERFACE JacobianBHSWME1D
    MODULE PROCEDURE JacobianBHSWME1D
END INTERFACE

INTERFACE JacobianSWLME1D
    MODULE PROCEDURE JacobianSWLME1D
END INTERFACE

INTERFACE JacobianMHSWME1D
    MODULE PROCEDURE JacobianMHSWME1D
END INTERFACE

INTERFACE JacobianPHSWME1D
    MODULE PROCEDURE JacobianPHSWME1D
END INTERFACE

INTERFACE JacobianPMHSWME1D
    MODULE PROCEDURE JacobianPMHSWME1D
END INTERFACE

INTERFACE JacobianInclinedHSWME1D
    MODULE PROCEDURE JacobianInclinedHSWME1D
END INTERFACE

INTERFACE JacobianReducedModelI
    MODULE PROCEDURE JacobianReducedModelI
END INTERFACE


INTERFACE QmatrixSWME1D
    MODULE PROCEDURE QmatrixSWME1D
END INTERFACE

INTERFACE QmatrixSWLME1D
    MODULE PROCEDURE QmatrixSWLME1D
END INTERFACE

INTERFACE ZeroMatrix
    MODULE PROCEDURE ZeroMatrix
END INTERFACE

INTERFACE ZeroSource
    MODULE PROCEDURE ZeroSource
END INTERFACE

INTERFACE SourceFuncSWME1D
    MODULE PROCEDURE SourceFuncSWME1D
END INTERFACE

INTERFACE SourceFuncSWE2D
    MODULE PROCEDURE SourceFuncSWE2D
END INTERFACE

INTERFACE SourceNewtonianSlip
    MODULE PROCEDURE SourceNewtonianSlip
END INTERFACE

INTERFACE SourceNewtonianManning
    MODULE PROCEDURE SourceNewtonianManning
END INTERFACE

INTERFACE SourceCoulombA
    MODULE PROCEDURE SourceCoulombA
END INTERFACE

INTERFACE SourceCoulombB
    MODULE PROCEDURE SourceCoulombB
END INTERFACE

INTERFACE SourceGranular
    MODULE PROCEDURE SourceGranular
END INTERFACE

INTERFACE SourceSavageHutter
    MODULE PROCEDURE SourceSavageHutter
END INTERFACE



INTERFACE SourceReducedModelI
    MODULE PROCEDURE SourceReducedModelI
END INTERFACE

INTERFACE FluxSWE1D
    MODULE PROCEDURE FluxSWE1D
END INTERFACE

INTERFACE FluxSWME1D
    MODULE PROCEDURE FluxSWME1D
END INTERFACE

INTERFACE FluxSWLME1D
    MODULE PROCEDURE FluxSWLME1D
END INTERFACE
!---------------------------------------
! two-dimensional
INTERFACE FluxSWME2DX
    MODULE PROCEDURE FluxSWME2DX
END INTERFACE

INTERFACE FluxSWME2DY
    MODULE PROCEDURE FluxSWME2DY
END INTERFACE

!-------------------------------------------------------------------------------!
! Note that what we call "Jacobian" are actually just the matrix A of the system:
! u_t + A(u) u_x  = 0;  
PUBLIC :: JacobianSWME1D
PUBLIC :: JacobianHSWME1D
PUBLIC :: JacobianBHSWME1D
PUBLIC :: JacobianSWLME1D
PUBLIC :: JacobianMHSWME1D
PUBLIC :: JacobianPHSWME1D
PUBLIC :: JacobianPMHSWME1D

PUBLIC :: JacobianInclinedHSWME1D
PUBLIC :: JacobianReducedModelI

PUBLIC :: QmatrixSWME1D
PUBLIC :: QmatrixSWLME1D
PUBLIC :: ZeroMatrix
PUBLIC :: ZeroSource
! 
PUBLIC :: SourceFuncSWME1D
PUBLIC :: SourceFuncSWE2D
PUBLIC :: SourceNewtonianSlip
PUBLIC :: SourceNewtonianManning
PUBLIC :: SourceCoulombA
PUBLIC :: SourceCoulombB
PUBLIC :: SourceGranular
PUBLIC :: SourceSavageHutter


PUBLIC :: SourceReducedModelI
PUBLIC :: FluxSWE1D
PUBLIC :: FluxSWME1D
PUBLIC :: FluxSWLME1D
!-------------------------------------------------------------------------------!
PUBLIC :: FluxSWME2DX
PUBLIC :: FluxSWME2DY
!-------------------------------------------------------------------------------!
! Pointers
PUBLIC :: SystemMatrix
PUBLIC :: RemainMatrix
PUBLIC :: Flux1D
PUBLIC :: Flux2DX
PUBLIC :: Flux2DY
PUBLIC :: SourceTerm
PUBLIC :: SourceSplit
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!


!-------------------------------------------------------------------------------!
! SHALLOW WATER MOMENT MODEL:
! Jacobian
!===============================================================================!
SUBROUTINE JacobianSWME1D(Prim, MatA)  ! One-dimensional
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms, Gravity, SubIdentity, VecDj
USE MOD_FiniteVolume_vars,ONLY: Acoef
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: MatA(1:nVar, 1:nVar)
!-------------------------------------------------------------------------------!


INTEGER  :: ii, jj, kk
REAL     :: h, um, alpha(1:nMoms), Psi

MatA = 0.0
h     = Prim(1)
um    = Prim(2)
alpha = Prim(3:nVar)

Psi = SUM(VecDj*alpha*alpha)
#ifndef TWODIM
!--------------------------------------------------------------------------------
MatA(1,2)      = 1.0
!--------------------------------------------------------------------------------
MatA(2,1)      = Gravity*h - um*um - Psi
MatA(2,2)      = 2.0 * um
MatA(2,3:nVar) = 2.0 * VecDj * alpha
!--------------------------------------------------------------------------------
MatA(3:nVar,1) = -2 * um * alpha
MatA(3:nVar,2) = 2.0 * alpha

DO ii = 1, nMoms
    DO jj = 1, nMoms
        ! Submatrix N x N              
        DO kk = 1,nMoms
             MatA(2 + ii, 1) = MatA(2 + ii,1) - Acoef(ii,jj,kk) * alpha(jj) * alpha(kk) 
             !--------------------------------------------------------------------------            
             MatA(2 + ii, 2 + jj) = MatA(2 + ii, 2+jj) + 2.0*Acoef(ii,jj,kk) * alpha(kk)
        END DO
    END DO
END DO
MatA(3:nVar,3:nVar) = MatA(3:nVar, 3:nVar) + 2.0 * um * SubIdentity   

#endif

END SUBROUTINE JacobianSWME1D
!===============================================================================!


!-------------------------------------------------------------------------------!
! Jacobian Hyperbolic SWME: HSWME
!===============================================================================!
SUBROUTINE JacobianHSWME1D(Prim, MatA)  ! One-dimensional
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms, Gravity, SubIdentity, VecDj
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: MatA(1:nVar, 1:nVar)
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj, kk
REAL             :: h, um, alpha(1:nMoms)

MatA = 0.0
h     = Prim(1)
um    = Prim(2)
alpha = Prim(3:nVar)

!-------------------------------------------------------------------------------!
! Acts only in 1D case:
#ifndef TWODIM
MatA(1,2) = 1.

MatA(2,1) = Gravity*h - um*um - (1./3.)*alpha(1)*alpha(1)
MatA(2,2) = 2.*um               
MatA(2,3) = (2./3.)*alpha(1)           

MatA(3,1) = -2.0 * um * alpha(1)
MatA(3,2) =  2.0 * alpha(1)

DO ii = 1, nMoms - 1
    ! Submatrix N x N
    MatA(ii + 2, ii + 3) = REAL(ii+2) * VecDj(ii+1) * alpha(1) ! upper-diagonal
    MatA(ii + 3, ii + 2) = REAL(ii)   * VecDj(ii)   * alpha(1) ! lower-diagonal

    MatA(4,1) = - (2.0/3.0) * alpha(1) * alpha(1)
END DO

MatA(3:nVar,3:nVar) = MatA(3:nVar, 3:nVar) + um * SubIdentity ! main lower block   
!-------------------------------------------------------------------------------!

#endif

END SUBROUTINE JacobianHSWME1D
!===============================================================================!


!-------------------------------------------------------------------------------!
! Jacobian beta Hyperbolic SWME: beta-HSWME
!===============================================================================!
SUBROUTINE JacobianBHSWME1D(Prim, MatA)  ! One-dimensional
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, VectorBeta
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: MatA(1:nVar, 1:nVar)
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj, kk

MatA = 0.0

#ifndef TWODIM
CALL JacobianHSWME1D(Prim, MatA)
MatA(nVar,:) = MatA(nVar,:) + VectorBeta  
#endif

END SUBROUTINE JacobianBHSWME1D
!===============================================================================!


!-------------------------------------------------------------------------------!
! Jacobian Linearized SWME: SWLME
!===============================================================================!
SUBROUTINE JacobianSWLME1D(Prim, MatA)  ! One-dimensional
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, Gravity, SubIdentity, VecDj, nMoms
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: MatA(1:nVar, 1:nVar)
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj, kk
REAL             :: h, um, alpha(1:nMoms), Psi

MatA  = 0.0
h     = Prim(1)
um    = Prim(2)
alpha = Prim(3:nVar)

Psi   = SUM(VecDj*alpha*alpha)
#ifndef TWODIM
MatA(1,2) = 1.                                 
MatA(3:nVar,1) = -2.* um * alpha;
MatA(3:nVar,2) =  2.* alpha;        
MatA(2,2)      =  2.* um;        

MatA(2,1)      = Gravity*h - um * um - Psi;
MatA(2,3:nVar) = 2. * VecDj * alpha;
! Submatrix N x N
MatA(3:nVar,3:nVar) = MatA(3:nVar, 3:nVar) + um * SubIdentity  
#endif

END SUBROUTINE JacobianSWLME1D
!===============================================================================!


!-------------------------------------------------------------------------------!
! Jacobian Moment regularization HSWME: MHSWME
!===============================================================================!
SUBROUTINE JacobianMHSWME1D(Prim, MatA)  ! One-dimensional
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, Gravity, SubIdentity, VecDj, nMoms
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: MatA(1:nVar, 1:nVar)
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj, kk
REAL             :: h, um, alpha(1:nMoms), Psi

MatA  = 0.0
h     = Prim(1)
um    = Prim(2)
alpha = Prim(3:nVar)

Psi   = SUM(VecDj*alpha*alpha)

#ifndef TWODIM
MatA(1,2)      = 1.

MatA(2,1)      = Gravity*h - um*um - Psi                               
MatA(2,2)      = 2.* um;  
MatA(2,3:nVar) = 2.* VecDj * alpha;

MatA(3,1)      = - 2.0* um * alpha(1)
MatA(3,2)      =   2.0* alpha(1)     
 
DO ii = 1, nMoms - 1
    ! Submatrix N x N
    MatA(ii + 2, ii + 3) = REAL(ii + 2) * VecDj(ii+1) * alpha(1) ! upper-diagonal
    MatA(ii + 3, ii + 2) = REAL(ii)     * VecDj(ii)   * alpha(1)          ! lower-diagonal
    MatA(4,1)            = - (2.0/3.0) * alpha(1) * alpha(1)
END DO
MatA(3:nVar,3:nVar) = MatA(3:nVar, 3:nVar) + um * SubIdentity ! main lower block   

#endif

END SUBROUTINE JacobianMHSWME1D
!===============================================================================!


!-------------------------------------------------------------------------------!
! Jacobian Primitive regularization HSWME: PHSWME
!===============================================================================!
SUBROUTINE JacobianPHSWME1D(Prim, MatA)  ! One-dimensional
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, Gravity, SubIdentity, VecDj, nMoms
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: MatA(1:nVar, 1:nVar)
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj, kk
REAL             :: h, um, alpha(1:nMoms)

MatA  = 0.0
h     = Prim(1)
um    = Prim(2)
alpha = Prim(3:nVar)

#ifndef TWODIM
MatA(1,2) = 1.

MatA(2,1) = Gravity*h - um*um - alpha(1)*alpha(1)/3.0;                                                             
MatA(2,2) = 2.* um;  
MatA(2,3) = (2.0/3.0)*alpha(1);

MatA(3,2)      = 2*alpha(1);
MatA(4:nVar,2) = alpha(2:nMoms);

MatA(3,1) = -2.0 * um * alpha(1) 
MatA(4,1) = -(1.0/3.0) * alpha(1) * alpha(1) 
 
DO ii = 1, nMoms - 1
    ! First column
    MatA(ii + 2,1) = MatA(ii + 2,1) - (REAL(ii+2)/REAL(2*ii+3)) * alpha(1) * alpha(ii+1)
    MatA(ii + 3,1) = MatA(ii + 3,1) -   (REAL(ii)/REAL(2*ii+1)) * alpha(1) * alpha(ii) - um * alpha(ii+1)
    
    ! Submatrix N x N
    MatA(ii + 2, ii + 3) = REAL(ii+2) * VecDj(ii+1) * alpha(1) ! upper-diagonal
    MatA(ii + 3, ii + 2) = REAL(ii)   * VecDj(ii)   * alpha(1) ! lower-diagonal
    
END DO
MatA(3:nVar,3:nVar) = MatA(3:nVar, 3:nVar) + um * SubIdentity ! main lower block   

#endif

END SUBROUTINE JacobianPHSWME1D
!===============================================================================!



!-------------------------------------------------------------------------------!
! Jacobian Primitive Moment regularization HSWME: PHSWME
!===============================================================================!
SUBROUTINE JacobianPMHSWME1D(Prim, MatA)  ! One-dimensional
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, Gravity, SubIdentity, VecDj, nMoms
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: MatA(1:nVar, 1:nVar)
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj, kk
REAL             :: h, um, alpha(1:nMoms), Psi

MatA  = 0.0
h     = Prim(1)
um    = Prim(2)
alpha = Prim(3:nVar)

Psi   = SUM(VecDj*alpha*alpha)
#ifndef TWODIM
MatA(1,2) = 1.

MatA(2,1) = Gravity*h - um*um - Psi                
MatA(2,2) = 2. * um;  
MatA(2,3:nVar) = 2. * VecDj * alpha;

MatA(3,2)      = 2. * alpha(1);
MatA(3,1) = -2.0 * um * alpha(1) 
 
DO ii = 1, nMoms - 1
    ! First column
    MatA(ii + 2,1) = MatA(ii + 2,1) - (REAL(ii+2)/REAL(2*ii+3)) * alpha(1) * alpha(ii+1)
    MatA(ii + 3,1) = MatA(ii + 3,1) -   (REAL(ii)/REAL(2*ii+1)) * alpha(1) * alpha(ii)   - um * alpha(ii+1)
    
    ! Submatrix N x N
    MatA(ii + 2, ii + 3) = REAL(ii+2) * VecDj(ii+1) * alpha(1) ! upper-diagonal
    MatA(ii + 3, ii + 2) = REAL(ii)   * VecDj(ii)   * alpha(1) ! lower-diagonal

    !--------------------------------------------
    MatA(4,1) = -(1.0/3.0) * alpha(1) * alpha(1) 
    MatA(4:nVar,2) = alpha(2:nMoms);   
    !--------------------------------------------     
END DO
MatA(3:nVar,3:nVar) = MatA(3:nVar, 3:nVar) + um * SubIdentity ! main lower block   

#endif

END SUBROUTINE JacobianPMHSWME1D
!===============================================================================!



!-------------------------------------------------------------------------------!
! Jacobian Inclined HSWME
!===============================================================================!
SUBROUTINE JacobianInclinedHSWME1D(Prim, MatA)  ! One-dimensional
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms, Gravity, SubIdentity, VecDj
USE MOD_FiniteVolume_vars,ONLY: AngleIncl, EpsilonRatio
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: MatA(1:nVar, 1:nVar)
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj, kk
REAL             :: h, um, alpha(1:nMoms)

MatA  = 0.0
h     = Prim(1)
um    = Prim(2)
alpha = Prim(3:nVar)

#ifndef TWODIM
MatA(1,2) = 1.

MatA(2,1) = EpsilonRatio*COS(AngleIncl)*h - um*um - (1./3.)*alpha(1)*alpha(1)                                                
MatA(2,2) = 2.* um                      
MatA(2,3) = (2./3.) * alpha(1)                  
MatA(3,1) = -2.0 * um * alpha(1)
MatA(3,2) =  2.0 * alpha(1)     
DO ii = 1, nMoms - 1
    ! Submatrix N x N   
    MatA(ii + 2, ii + 3) = REAL(ii+2) * VecDj(ii+1) * alpha(1); ! upper-diagonal
    MatA(ii + 3, ii + 2) = REAL(ii)   * VecDj(ii)   * alpha(1); ! lower-diagonal
    
    MatA(4,1) = -(2.0/3.0) * alpha(1) * alpha(1)
END DO
MatA(3:nVar,3:nVar) = MatA(3:nVar, 3:nVar) + um * SubIdentity ! main   
!-------------------------------------------------------------------------------!
#endif

END SUBROUTINE JacobianInclinedHSWME1D
!===============================================================================!


!-------------------------------------------------------------------------------!
! Jacobian Reduced Model I
!===============================================================================!
SUBROUTINE JacobianReducedModelI(Prim, MatA)  ! One-dimensional
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, Gravity, nMoms, nMomsOut
USE MOD_FiniteVolume_vars,ONLY: ParameterLambda
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: MatA(1:nVar, 1:nVar)
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj, kk

MatA = 0.0

#ifndef TWODIM
MatA(1,2)  = 1. 

IF (nMomsOut .EQ. 0) THEN

    MatA(2,1) = - Prim(2) * Prim(2) + Gravity * Prim(1) 
    MatA(2,2) = 2.*Prim(2)

ELSE IF (nMomsOut .EQ. 1) THEN

    MatA(2,1) = Prim(2) * Prim(2) * ( -1.0  +  Prim(1) * Prim(1)/(48.0*ParameterLambda*ParameterLambda))     &
                + Gravity * Prim(1) - 0.*(Gravity * Prim(1) * Prim(1) * Prim(1)/(48. * ParameterLambda * ParameterLambda))                            
    MatA(2,2) = 2.*Prim(2)  +  Prim(2) * Prim(1) * Prim(1)/( 24. * ParameterLambda * ParameterLambda)
                
                
ELSE IF ((nMomsOut .GT. 1) .AND. (nMomsOut .LT. 8)) THEN

    MatA(2,1) = Prim(2) * Prim(2) * ( -1.0  +  Prim(1) * Prim(1)/(45.0*ParameterLambda*ParameterLambda))     &
                + Gravity * Prim(1) - 0.*(Gravity * Prim(1) * Prim(1) * Prim(1)/(45. * ParameterLambda * ParameterLambda))
    MatA(2,2) = 2.*Prim(2)*(1.0 + Prim(1) * Prim(1)/( 45. * ParameterLambda * ParameterLambda))
    
END IF
!only works up to N = 7                         
#endif

END SUBROUTINE JacobianReducedModelI
!===============================================================================!


!-------------------------------------------------------------------------------!
! Non local matrices Q 
!===============================================================================!
SUBROUTINE QmatrixSWME1D(Prim, MatB)  ! One-dimensional
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms, SubIdentity, Bcoef 
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: MatB(1:nVar, 1:nVar)
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj, kk
REAL             :: h, um, alpha(1:nMoms)

MatB  = 0.0
h     = Prim(1)
um    = Prim(2)
alpha = Prim(3:nVar)

#ifndef TWODIM
DO ii = 1, nMoms
    DO jj = 1, nMoms
        MatB(2 + ii, 2 + jj) = SUM(Bcoef(ii,jj,1:nMoms) * alpha);
    END DO
END DO
MatB(3:nVar,3:nVar) = um * SubIdentity - MatB(3:nVar,3:nVar) 
MatB = - MatB

#endif




END SUBROUTINE QmatrixSWME1D
!===============================================================================!


!-------------------------------------------------------------------------------!
! Non local matrix Q - Linear
!===============================================================================!
SUBROUTINE QmatrixSWLME1D(Prim, MatB)  ! One-dimensional
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, SubIdentity
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: MatB(1:nVar, 1:nVar)

MatB = 0.0
#ifndef TWODIM
MatB(3:nVar, 3:nVar) = Prim(2)*SubIdentity 
#endif
  
END SUBROUTINE QmatrixSWLME1D
!===============================================================================!


!-------------------------------------------------------------------------------!
! SOURCE TERMS:
!===============================================================================!
SUBROUTINE SourceFuncSWME1D(Prim, SourceSW)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms, ParameterNu, ParameterLambda
USE MOD_FiniteVolume_vars,ONLY: Ccoef, VecDj
!-------------------------------------------------------------------------------!
IMPLICIT NONE
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: SourceSW(1:nVar)
REAL    :: lmbda, nu
INTEGER :: ii, jj
REAL    :: SP(1:nVar), slip_length, PP(1:nMoms+1),  matc(1:nMoms+1,nMoms)

SourceSW = 0.0
#ifndef TWODIM
nu    = ParameterNu     ! nu
lmbda = ParameterLambda ! lambda
SourceSW(1) = 0.
SourceSW(2) = - (nu/lmbda)*(Prim(2) + SUM(Prim(3:nVar)))
DO ii = 1, nMoms                        
    SourceSW(2 + ii) = REAL(2*ii+1) * ( SourceSW(2)  - (nu/Prim(1)) * SUM(Ccoef(ii,1:nMoms)*Prim(3:nVar)))    
END DO
#endif

END SUBROUTINE SourceFuncSWME1D
!===============================================================================!


!===============================================================================!
SUBROUTINE SourceNewtonianSlip(Prim, SourceSW)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms, ParameterNu, ParameterLambda
USE MOD_FiniteVolume_vars,ONLY: Ccoef, VecDj
USE MOD_FiniteVolume_vars,ONLY: AngleIncl, EpsilonRatio, IndexPass, dBath, Gravity
!-------------------------------------------------------------------------------!
IMPLICIT NONE
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: SourceSW(1:nVar)
!-------------------------------------------------------------------------------!
REAL    :: nu, lmbda, Auxsum, TauXZB, Lconst, Hconst, Uconst, Rho
INTEGER :: ii

Lconst = 10.0
Hconst = EpsilonRatio*Lconst
Uconst = SQRT(Gravity*Lconst)
Rho    = 1200.0

SourceSW = 0.0

lmbda = ParameterLambda/Hconst  ! lambda dimensionless
nu    = ParameterNu*Uconst/(Rho*Gravity*COS(AngleIncl)*Hconst*Hconst) ! nu

#ifndef TWODIM
Auxsum = SUM(PRIM(2:nVar))
TauXZB = (nu/lmbda)*Auxsum
SourceSW(1) = 0.
SourceSW(2) = SIN(AngleIncl) * Prim(1) - COS(AngleIncl) * TauXZB &
               - EpsilonRatio * COS(AngleIncl) * Prim(1) * dBath(IndexPass(1), IndexPass(2))

DO ii = 1, nMoms
    SourceSW(2 + ii) = -REAL(2*ii+1) * COS(AngleIncl)*(TauXZB + (nu/Prim(1))*SUM(Ccoef(ii,1:nMoms)*Prim(3:nVar)) )
END DO
#endif
!-------------------------------------------------------------------------------!
END SUBROUTINE SourceNewtonianSlip
!===============================================================================!


!===============================================================================!
SUBROUTINE SourceNewtonianManning(Prim, SourceSW)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms, ParameterNu, ParameterLambda
USE MOD_FiniteVolume_vars,ONLY: Ccoef, VecDj, AuxiliaryParameter
USE MOD_FiniteVolume_vars,ONLY: AngleIncl, EpsilonRatio, IndexPass, dBath, Gravity
!-------------------------------------------------------------------------------!
IMPLICIT NONE
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: SourceSW(1:nVar)
!-------------------------------------------------------------------------------!
REAL    :: nu, lmbda, Auxsum, TauXZB, nConstant2, Lconst, Hconst, Uconst, Rho, nConstant0
INTEGER :: ii
!-------------------------------------------------------------------------------!
SourceSW = 0.0


Lconst = 10.0
Hconst = EpsilonRatio*Lconst
Uconst = SQRT(Gravity*Lconst)
Rho    = 1200.0

!nConstant0 = 0.0165 # this is the value from papers
nConstant0 = AuxiliaryParameter

!-----------------------------------------------------------------------
nu     = ParameterNu*Uconst/(Rho*Gravity*COS(AngleIncl)*Hconst*Hconst)   ! nu
! Power 2 of n constant
nConstant2 = nConstant0**2*Uconst*Uconst/(Hconst**(4.0/3.0)*COS(AngleIncl))  ! n^2 is the first number here

Auxsum = SUM(PRIM(2:nVar))

TauXZB = (nConstant2/((Prim(1))**(1./3.)))*ABS(Auxsum)*Auxsum

#ifndef TWODIM
SourceSW(1) = 0.
SourceSW(2) = SIN(AngleIncl) * Prim(1) - COS(AngleIncl) * TauXZB &
              - EpsilonRatio * COS(AngleIncl) * Prim(1) * dBath(IndexPass(1), IndexPass(2))
DO ii = 1, nMoms
    SourceSW(2 + ii) = -REAL(2*ii+1) * COS(AngleIncl)*(TauXZB + &
                        (nu/Prim(1))*DOT_PRODUCT(Ccoef(ii,1:nMoms),Prim(3:nVar)) )
END DO
#endif
!-------------------------------------------------------------------------------!
END SUBROUTINE SourceNewtonianManning
!===============================================================================!


!===============================================================================!
SUBROUTINE SourceSavageHutter(Prim, SourceSW)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms, ParameterNu, ParameterLambda, PI
USE MOD_FiniteVolume_vars,ONLY: Ccoef, VecDj, AuxiliaryParameter
USE MOD_FiniteVolume_vars,ONLY: AngleIncl, EpsilonRatio, IndexPass, dBath
USE MOD_FiniteVolume_vars,ONLY: QuadSourceWeights, QuadSourceNodes, nGSource
USE MOD_PhysicsFrame,ONLY: LegendrePoly, LegendrePolyNormalized01
!-------------------------------------------------------------------------------!
IMPLICIT NONE
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: SourceSW(1:nVar)
!-------------------------------------------------------------------------------!
REAL    :: nu, lmbda, Auxsum, TauXZB, TauI(1:nMoms), ANGLE1, ANGLE2, SGNM, SGNM0
REAL    :: AuxDphi(1:nMoms), AuxPhi(1:nMoms), Qvalue, PolyN, dPolyN, s0aux, Auxprod
INTEGER :: ii, iGG, jj
!-------------------------------------------------------------------------------!

SourceSW = 0.0
AuxDphi  = 0.0

ANGLE1 = AuxiliaryParameter*PI/180.0  ! friction angle
ANGLE2 = 20.0*PI/180.0                ! internal angle

SGNM0  = SIGN(1.0, SUM(Prim(2:nVar)))
IF ( SUM(Prim(2:nVar)) .EQ. 0.0 ) THEN
   SGNM0 = 0.0
ENDIF
TauXZB = Prim(1)*TAN(ANGLE1)*SGNM0
! Computation of the \mathcal{T}_i
TauI = 0.0
TauI = -Prim(1)*TAN(ANGLE2)

!print*, TauXZB - TauI(1)
#ifndef TWODIM
SourceSW(1) = 0.
SourceSW(2) = SIN(AngleIncl)*Prim(1) - COS(AngleIncl)*TauXZB &
              - EpsilonRatio* COS(AngleIncl) * Prim(1) * dBath(IndexPass(1), IndexPass(2))
!-------------------------------------------------------------------------------!
DO ii = 1, nMoms
    SourceSW(2 + ii) =  -REAL(2*ii+1) * COS(AngleIncl)*(TauXZB + TauI(ii))
END DO
#endif
!-------------------------------------------------------------------------------!
END SUBROUTINE SourceSavageHutter
!===============================================================================!


!===============================================================================!
SUBROUTINE SourceCoulombA(Prim, SourceSW)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms
USE MOD_FiniteVolume_vars,ONLY: PI, Gravity
USE MOD_FiniteVolume_vars,ONLY: Ccoef, VecDj, AuxiliaryParameter
USE MOD_FiniteVolume_vars,ONLY: AngleIncl, EpsilonRatio, IndexPass, dBath
USE MOD_FiniteVolume_vars,ONLY: QuadSourceWeights, QuadSourceNodes, nGSource
USE MOD_PhysicsFrame,ONLY: LegendrePoly, LegendrePolyNormalized01
!-------------------------------------------------------------------------------!
IMPLICIT NONE
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: SourceSW(1:nVar)
!-------------------------------------------------------------------------------!
REAL    :: CTE3, Auxsum, TauXZB, TauI(1:nMoms), ANGLE1, SGNM, SGNM0
REAL    :: Lconst, Hconst, Uconst, Rho, MU0, Auxprod
REAL    :: AuxDphi(1:nMoms), Qvalue, PolyN, dPolyN, s0aux
INTEGER :: ii, iGG, jj
!-------------------------------------------------------------------------------!
SourceSW = 0.0
AuxDphi  = 0.0

Lconst = 10.0
Hconst = EpsilonRatio*Lconst
Uconst = SQRT(Gravity*Lconst)
Rho    = 1200.0

! Be aware that an angle mu0 with 29 and delta with 15 gives a backwave
MU0    = TAN(20.0*PI/180.0)
ANGLE1 = AuxiliaryParameter*PI/180.0 ! friction angle
Auxsum = SUM(Prim(2:nVar))

! Relation between tan(theta) <= mu_s (Theorem 1, Garres2021b)
SGNM0  = SIGN(1.0, SUM(Prim(2:nVar)))
IF ( SUM(Prim(2:nVar)) .EQ. 0.0 ) THEN
   SGNM0 = 0.0
ENDIF
TauXZB = Prim(1)*TAN(ANGLE1)*SGNM0
! Computation of the \mathcal{T}_i
TauI = 0.0
TauI = -Prim(1)*MU0

!-------------------------------------------------------------------------------!
! (48) FORMULA
#ifndef TWODIM
SourceSW(1) = 0.
SourceSW(2) = SIN(AngleIncl)*Prim(1) - COS(AngleIncl)*TauXZB &
               - EpsilonRatio* COS(AngleIncl) * Prim(1) * dBath(IndexPass(1), IndexPass(2))
!-------------------------------------------------------------------------------!
DO ii = 1, nMoms
    SourceSW(2 + ii) = - REAL(2*ii+1) * COS(AngleIncl)*(TauXZB + TauI(ii)) 
END DO
#endif
!-------------------------------------------------------------------------------!
END SUBROUTINE SourceCoulombA
!===============================================================================!






!===============================================================================!
SUBROUTINE SourceCoulombB(Prim, SourceSW)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms
USE MOD_FiniteVolume_vars,ONLY: PI, Gravity
USE MOD_FiniteVolume_vars,ONLY: Ccoef, VecDj, AuxiliaryParameter
USE MOD_FiniteVolume_vars,ONLY: AngleIncl, EpsilonRatio, IndexPass, dBath
USE MOD_FiniteVolume_vars,ONLY: QuadSourceWeights, QuadSourceNodes
USE MOD_PhysicsFrame,ONLY: LegendrePoly, LegendrePolyNormalized01
!-------------------------------------------------------------------------------!
IMPLICIT NONE
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: SourceSW(1:nVar)
!-------------------------------------------------------------------------------!
SourceSW = 0.0
!-------------------------------------------------------------------------------!
END SUBROUTINE SourceCoulombB
!===============================================================================!




!===============================================================================!
SUBROUTINE SourceGranular(Prim, SourceSW)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms, ParameterNu, ParameterLambda
USE MOD_FiniteVolume_vars,ONLY: ParameterViscosity, Gravity, PI
USE MOD_FiniteVolume_vars,ONLY: Ccoef, VecDj, AuxiliaryParameter
USE MOD_FiniteVolume_vars,ONLY: AngleIncl, EpsilonRatio, IndexPass, dBath, nGSource
USE MOD_FiniteVolume_vars,ONLY: QuadSourceWeights, QuadSourceNodes, LegenPoly, DxLegenPoly
USE MOD_PhysicsFrame,ONLY: LegendrePolyNormalized01
!-------------------------------------------------------------------------------!
IMPLICIT NONE
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: SourceSW(1:nVar)
!-------------------------------------------------------------------------------!
REAL    :: nu, lmbda, CTE3, TauXZB, MUS, MU2, RHO, RHOS, VAL0, VAL1, TauI(1:nMoms), I0const
REAL    :: Lconst, Hconst, Uconst, dsconst, ANGLE0
REAL    :: Auxprod, AuxDphi(1:nMoms), Qvalue, SGNM, SGNM0, PolyN, dPolyN, Auxprod0
REAL    :: s0aux, Auxsum
INTEGER :: ii, iGG, jj
!-------------------------------------------------------------------------------!

SourceSW = 0.0; AUXPROD  = 0.0; AuxDphi  = 0.0

Lconst = 10.0
Hconst = EpsilonRatio*Lconst
Uconst = SQRT(Gravity*Lconst)
Rho    = 1200.0

nu     = ParameterNu        ! nu
lmbda  = ParameterLambda    ! lambda

I0const = 0.279
dsconst = 0.7

!mu = 0.48 + (0.73 - 0.48)/(0.279 + I)*I
MUS  = 0.48
MU2  = 0.73
!RHO  = 1000.0 
!RHOS = 1580.0 
RHO  = 1550.0 
RHOS = 2500.0 
VAL0 = (dsconst/ Hconst) * SQRT(RHOS/(RHO*EpsilonRatio*COS(AngleIncl)))  ! c_I constant
VAL1 = (I0const*PRIM(1)**(3/2))/VAL0                                     ! Constant C in the manuscript
!---------------------------------------------------------------------------------------

ANGLE0 = 20.0*PI/180.0    ! friction angle
SGNM0  = SIGN(1.0,SUM(Prim(2:nVar)))
IF (SUM(Prim(2:nVar)) .EQ. 0.0) THEN
   SGNM0 = 0.0
ENDIF
TauXZB = Prim(1)*TAN(ANGLE0)*SGNM0


lmbda  = ParameterLambda/Hconst  ! lambda dimensionless
nu     = ParameterNu*Uconst/(Rho*Gravity*COS(AngleIncl)*Hconst*Hconst) ! nu
Auxsum = SUM(PRIM(2:nVar))
TauXZB = (nu/lmbda)*Auxsum

TauI = 0.0
! Quadrature rule:
DO iGG = 1,nGSource
    s0aux = QuadSourceNodes(iGG)
    !------------------------------------------------------------
    Auxprod = SUM(Prim(3:nVar)*DxLegenPoly(iGG,:))   
    SGNM = SIGN(1.0,Auxprod)
    IF (Auxprod .EQ. 0.0) THEN
        SGNM = 0.0
    END IF
    Qvalue = (1.0 - s0aux)*(MUS + ((MU2 - MUS)*ABS(Auxprod))/(ABS(Auxprod) + VAL1*SQRT(1.0 - s0aux))) * SGNM    
    TauI = TauI + QuadSourceWeights(iGG) * Prim(1) * Qvalue * DxLegenPoly(iGG,:)
END DO
!-------------------------------------------------------------------------------!
#ifndef TWODIM
SourceSW(1) = 0.
SourceSW(2) = SIN(AngleIncl)*Prim(1) - COS(AngleIncl)*TauXZB &
              - EpsilonRatio* COS(AngleIncl) * Prim(1) * dBath(IndexPass(1), IndexPass(2))
!-------------------------------------------------------------------------------!
DO ii = 1, nMoms
    SourceSW(2 + ii) = - REAL(2*ii+1) * COS(AngleIncl)*(TauXZB + TauI(ii))
END DO
#endif
!-------------------------------------------------------------------------------!
END SUBROUTINE SourceGranular
!===============================================================================!




















!===============================================================================!
SUBROUTINE SourceFuncSWE2D(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars, ONLY: U,V,SSFriction
USE MOD_FiniteVolume_vars, ONLY: nVar, nXs, nYs, nDims, nGPsX, nGPsY
USE MOD_FiniteVolume_vars, ONLY: MeshBary, MeshGP, WeightsGP, MESH_DX
USE MOD_FiniteVolume_vars, ONLY: nGhostsX, nGhostsY 
USE MOD_FiniteVolume_vars, ONLY: ReconstrFlag, BathymetryFlag
USE MOD_PhysicsFrame,      ONLY: Bathymetry_X, Bathymetry_Y
USE MOD_Reconstruction,    ONLY: MUSCL, WENO3_SecondSweep, WENO5_SecondSweep
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: t
!-------------------------------------------------------------------------------!
REAL               :: SW(1:nVar, 1:nGPsX, 1:nGPsY, 1:nXs, 1:nYs)
REAL               :: Vtemp(1:nVar, -nGhostsX:nXs + nGhostsX + 1, 1:nYs, 1:nGPsX)
REAL               :: Vtemp2(1:nVar, 1:nGPsX, 1:nGPsY, 1:nXs, 1:nYs)
INTEGER            :: ii, jj, iVar, iGP, jGP
CHARACTER(LEN=255) :: ErrorMessage

SSFriction = 0.0 ! S(1:nVar,nXs,nYs)
!int_{iixjj} S(1:nVar,ii,jj) dxdy

IF (BathymetryFlag .GT. 0) THEN
    SELECT CASE (ReconstrFlag)
    CASE (1, 2)
        DO jj = 1, nYs
            DO ii = 1, nXs
                SW(:, nGPsX, nGPsY, ii, jj) = SourceFunc(U(:, ii, jj), MeshBary(:, ii, jj))
            END DO
        END DO       
    !----------------------------------------------------
    CASE (3)
    
        DO iVar = 1, nVar
            DO jj = 1, nYs
                DO ii = -nGhostsX, nXs + nGhostsX + 1
                    CALL WENO3_SecondSweep(U(iVar, ii, jj - nGhostsY:jj + nGhostsY),  Vtemp(iVar, ii, jj,:))
                END DO
            END DO
        END DO
        
        DO jj = 1, nYs
            DO ii = 1, nXs
                DO iGP = 1, nGPsX
                    DO iVar = 1, nVar
                        CALL WENO3_SecondSweep(Vtemp(iVar, ii - nGhostsX:ii + nGhostsX, jj, iGP), Vtemp2(iVar,:, iGP, ii, jj))
                    END DO
                    DO jGP = 1, nGPsY
                        SW(:, jGP, iGP, ii, jj) = SourceFunc(Vtemp2(:, jGP, iGP, ii, jj), MeshGP(:, ii, jj, jGP, iGP))
                    END DO
                END DO
            END DO
        END DO

    !----------------------------------------------------        
    CASE (4)
    
        DO iVar = 1, nVar
            DO jj = 1, nYs
                DO ii = -nGhostsX,  nXs + nGhostsX + 1
                    CALL WENO5_SecondSweep(U(iVar, ii, jj - nGhostsY:jj + nGhostsY), Vtemp(iVar, ii, jj,:))
                END DO
            END DO
        END DO
        
        DO jj = 1, nYs
            DO ii = 1, nXs
                DO iGP = 1, nGPsX
                    DO iVar = 1, nVar
                        CALL WENO5_SecondSweep(Vtemp(iVar, ii - nGhostsX:ii + nGhostsX, jj, iGP), Vtemp2(iVar,:, iGP, ii, jj))
                    END DO
                    
                    DO jGP = 1, nGPsY
                        SW(:, jGP, iGP, ii, jj) = SourceFunc(Vtemp2(:, jGP, iGP, ii, jj), MeshGP(:, ii, jj, jGP, iGP))
                    END DO
                END DO
            END DO
        END DO

    !----------------------------------------------------
    CASE DEFAULT
        ErrorMessage = "Reconstruction not implemented"
        WRITE (*, *) ErrorMessage
        STOP
    END SELECT

    DO jj = 1, nYs
        DO ii = 1, nXs
            DO iGP = 1, nGPsX
                DO jGP = 1, nGPsY
                    SSFriction(:, ii, jj) = SSFriction(:, ii, jj) + WeightsGP(iGP, jGP)*SW(:, iGP, jGP, ii, jj)
                END DO
            END DO
        END DO
    END DO

END IF

CONTAINS

    FUNCTION SourceFunc(Q,X) RESULT(Source_SW)
    !-------------------------------------------------------------------------------!
    USE MOD_FiniteVolume_vars,ONLY: nDims, nVar, Gravity
    IMPLICIT NONE
    REAL, DIMENSION(1:nVar),  INTENT(IN)  :: Q 
    REAL, DIMENSION(1:nDims), INTENT(IN)  :: X
    REAL, DIMENSION(1:nVar) :: Source_SW 
    !-------------------------------------------------------------------------------!
    INTEGER :: ii
    !-------------------------------------------------------------------------------!
    Source_SW(1) = 0.
    Source_SW(2) = -Gravity*Q(1)*Bathymetry_X(X) 
#ifdef TWODIM
    Source_SW(3) = -Gravity*Q(1)*Bathymetry_Y(X)
#endif
    !-------------------------------------------------------------------------------!
    END FUNCTION SourceFunc


!-------------------------------------------------------------------------------!
END SUBROUTINE SourceFuncSWE2D
!===============================================================================!


!===============================================================================!
SUBROUTINE SourceReducedModelI(Prim, SourceSW)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms, ParameterNu, ParameterLambda, nMomsOut
!-------------------------------------------------------------------------------!
IMPLICIT NONE
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: SourceSW(1:nVar)
!-------------------------------------------------------------------------------!
REAL    :: nu
REAL    :: lmbda
INTEGER :: ii

SourceSW = 0.0

! Acts only in 1D case:
#ifndef TWODIM
!-------------------------------------------------------------------------------!
! Here Prim is considered as the primitive (h, u, alpha1, dots, alphaN)
nu = ParameterNu     ! nu
lmbda = ParameterLambda ! lambda
 
SourceSW(1) = 0.

IF (nMomsOut .EQ. 0) THEN
    SourceSW(2) = - nu*Prim(2)/lmbda
ELSE IF (nMomsOut .EQ. 1) THEN
    SourceSW(2) = - nu*Prim(2)/lmbda * (1.0 - Prim(1)/(4.0*lmbda) + 0.*Prim(1)*Prim(1)/(24.0*lmbda*lmbda))
ELSE IF ((nMomsOut .GT. 1) .AND. (nMomsOut .LT. 8) ) THEN 
    SourceSW(2) = - nu*Prim(2)/lmbda * (1.0 - Prim(1)/(3.0*lmbda) + 0.*4.0*Prim(1)*Prim(1)/(45.0*lmbda*lmbda))
ELSE 
    SourceSW(2) = - nu*Prim(2)/lmbda * (1.0 - Prim(1)/(3.0**lmbda) + 0.*4.0*Prim(1)*Prim(1)/(45.0*lmbda*lmbda))
ENDIF

!This function only works for values between nMoms<8.
!-------------------------------------------------------------------------------!
#endif
!-------------------------------------------------------------------------------!
END SUBROUTINE SourceReducedModelI
!===============================================================================!


!===============================================================================!
SUBROUTINE FluxSWE1D(Prim,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, Gravity, Hmin, MIN_SPEED
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
REAL             :: h, vx, vy

#ifdef TWODIM
h  = Prim(1)
vx = Prim(2)
vy = Prim(3)

Flux(1) = h*vx
Flux(2) = h*vx**2 + 0.5*Gravity*h**2
Flux(3) = h*vx*vy
#endif
!-------------------------------------------------------------------------------!
END SUBROUTINE FluxSWE1D
!===============================================================================!


!-------------------------------------------------------------------------------!
! FLUXES:
!===============================================================================!
SUBROUTINE FluxSWME1D(Prim,Flux) ! Only Conservative Flux
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms, Gravity, VecDj, Acoef
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
REAL     :: h
INTEGER  :: ii, jj, kk
!-------------------------------------------------------------------------------!
Flux = 0.0
!-------------------------------------------------------------------------------!
! Acts only in 1D case:
#ifndef TWODIM

h       = Prim(1)
Flux(1) = h*Prim(2)
Flux(2) = h*(Prim(2))**2 + 0.5*Gravity*h**2 + h*DOT_PRODUCT(VecDj, Prim(3:nVar) * Prim(3:nVar))
DO ii = 1, nMoms
    Flux(2+ii) = 2.*h*Prim(2)*Prim(2+ii)
    DO jj = 1, nMoms
        DO kk = 1, nMoms
            Flux(2+ii) = Flux(2+ii) + h * Acoef(ii,jj,kk) * Prim(2+jj) * Prim(2+kk)
        END DO
    END DO
END DO

#endif

END SUBROUTINE FluxSWME1D
!===============================================================================!


!===============================================================================!
SUBROUTINE FluxSWLME1D(Prim,Flux)   ! Only Conservative Flux
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, Gravity, VecDj
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
REAL     :: h
INTEGER  :: ii, jj, kk

Flux = 0.0

! Acts only in 1D case:
#ifndef TWODIM

h       = Prim(1)
Flux(1) = h*Prim(2)
Flux(2) = h*(Prim(2))**2 + 0.5*Gravity*h**2 + h*DOT_PRODUCT(VecDj, Prim(3:nVar) * Prim(3:nVar))
Flux(3:nVar) = 2.*h*Prim(2)*Prim(3:nVar)

#endif

END SUBROUTINE FluxSWLME1D
!===============================================================================!


!===============================================================================!
SUBROUTINE FluxSWME2DX(Prim,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms, Gravity, Acoef, VecDj
USE MOD_FiniteVolume_vars,ONLY: Hmin, MIN_SPEED
IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
REAL             :: h, vx, vy
INTEGER          :: jj, ii, kk

Flux = 0.0
#ifdef TWODIM
h  = Prim(1)
vx = Prim(2)
vy = Prim(3)

Flux(1) = h*vx
Flux(2) = h*vx**2 + 0.5*Gravity*h**2
Flux(3) = h*vx*vy

! alphas and betas 
DO ii = 1, nMoms ! j
    Flux(2) = Flux(2) + (1./(2.*(REAL(ii)) + 1.))*h*Prim( 2*ii-1 + 3)**2
    Flux(3) = Flux(3) + (1./(2.*(REAL(ii)) + 1.))*h*Prim( 2*ii-1 + 3)*Prim( 2*ii + 3 )
    
    Flux(2*ii-1  + 3) = h*vx*Prim(2*ii-1 + 3) + h*vx*Prim(2*ii-1 + 3) 
    Flux(2*ii    + 3) = h*vy*Prim(2*ii-1 + 3) + h*vx*Prim(2*ii   + 3)
    
    DO jj = 1,nMoms ! i
        DO kk = 1, nMoms ! k
             Flux(2*ii-1 + 3) = Flux(2*ii-1 + 3) + Acoef(ii,jj,kk) * h*Prim(2*jj - 1 + 3) * Prim(2*kk - 1 + 3)
             Flux(2*ii   + 3) = Flux(2*ii   + 3) + Acoef(ii,jj,kk) * h*Prim(2*jj + 3)     * Prim(2*kk - 1 + 3)
        END DO
    END DO
END DO
#endif

END SUBROUTINE FluxSWME2DX
!===============================================================================!


!===============================================================================!
SUBROUTINE FluxSWME2DY(Prim,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nMoms, Gravity, Acoef, VecDj
USE MOD_FiniteVolume_vars,ONLY: Hmin, MIN_SPEED
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
REAL             :: h, vx, vy
INTEGER          :: ii, jj, kk
!-------------------------------------------------------------------------------!
Flux = 0.0
#ifdef TWODIM
h  = Prim(1)
vx = Prim(2)
vy = Prim(3)

Flux(1) = h*vy
Flux(2) = h*vx*vy
Flux(3) = h*vy**2 + 0.5*Gravity*h**2
!-------------------------------------------------------------------------------!
!!! alphas and betas 
DO ii = 1, nMoms
    Flux(2) = Flux(2) + (1./(2.*(REAL(ii)) + 1.))*h*Prim( 2*ii-1 + 3)*Prim( 2*ii + 3 )
    Flux(3) = Flux(3) + (1./(2.*(REAL(ii)) + 1.))*h*Prim( 2*ii   + 3)**2
    
    Flux(2*ii-1  + 3) = h*vx*Prim(2*ii + 3) + h*vy*Prim(2*ii-1 + 3) 
    Flux(2*ii    + 3) = h*vy*Prim(2*ii + 3) + h*vy*Prim(2*ii   + 3)
    
    DO jj = 1,nMoms ! i
        DO kk = 1, nMoms ! k
             Flux(2*ii-1 + 3) = Flux(2*ii-1 + 3) + Acoef(ii,jj,kk) * h*Prim(2*jj - 1 + 3) * Prim(2*kk + 3) ! alpha * beta
             Flux(2*ii   + 3) = Flux(2*ii   + 3) + Acoef(ii,jj,kk) * h*Prim(2*jj + 3)     * Prim(2*kk + 3) ! beta  * beta
        END DO
    END DO
END DO
!-------------------------------------------------------------------------------!
#endif
END SUBROUTINE FluxSWME2DY
!===============================================================================!
! 
! 
! 
!  
!!===============================================================================!
!FUNCTION SourceFuncSWE(Q,X) RESULT(Source_SW)
!!-------------------------------------------------------------------------------!
!USE MOD_FiniteVolume_vars,ONLY: nDims, nVar, Gravity
!IMPLICIT NONE
!REAL, DIMENSION(1:nVar),  INTENT(IN)  :: Q 
!REAL, DIMENSION(1:nDims), INTENT(IN)  :: X
!REAL, DIMENSION(1:nVar) :: Source_SW 
!INTEGER                 :: jj, ii, kk
!REAL                    :: const, viscosity, rho
!!-------------------------------------------------------!
!rho       = 1000.
!viscosity = 0.001
!const     = viscosity/(rho*Q(1))
!!-------------------------------------------------------!

!Source_SW(1) = 0.
!Source_SW(2) = -Gravity*Q(1)*Bathymetry_X(X) 
!Source_SW(3) = -Gravity*Q(1)*Bathymetry_Y(X) 

!DO jj = 1, nVar

!    Source_SW(2*jj-1 + 3) = Q(2)*XXXXX
!    Source_SW(2*jj   + 3) = Q(3)*XXXXX


!    DO ii = 1,nMoms ! i
!        DO kk = 1, nMoms
!            Source_SW(2*jj-1 + 3) = Source_SW(2*jj-1 + 3) - Bcoef(ii,jj,kk) * Q(2*jj-1 + 3) * XXXXXX
!            Source_SW(2*jj   + 3) = Source_SW(2*jj   + 3) - Bcoef(ii,jj,kk) * Q(2*jj   + 3) * XXXXXX
!        END DO
!        
!        Source_SW(2*jj-1 + 3) = Source_SW(2*jj-1 + 3) + const * Ccoef(ii,jj) * Q(2*jj-1 + 3)
!        Source_SW(2*jj   + 3) = Source_SW(2*jj   + 3) + const * Ccoef(ii,jj) * Q(2*jj   + 3)
!        
!    END DO

!END DO
!! TODO: define Ccoef and Bcoef

!!-------------------------------------------------------------------------------!
!END FUNCTION SourceFuncSWE
!!===============================================================================!


SUBROUTINE ZeroMatrix(Prim, MatA)
    USE MOD_FiniteVolume_vars,ONLY: nVar
    REAL,INTENT(IN)  :: Prim(1:nVar)
    REAL,INTENT(OUT) :: MatA(1:nVar, 1:nVar)            
    MatA = 0.0
END SUBROUTINE ZeroMatrix


SUBROUTINE ZeroSource(Prim, SourceSW)
    USE MOD_FiniteVolume_vars,ONLY: nVar
    REAL,INTENT(IN)  :: Prim(1:nVar)
    REAL,INTENT(OUT) :: SourceSW(1:nVar)
    SourceSW = 0.0
END SUBROUTINE ZeroSource

    
!===============================================================================!
END MODULE MOD_MomentModels
!-------------------------------------------------------------------------------!
