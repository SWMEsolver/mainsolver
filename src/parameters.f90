!===============================================================================!
MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE InitializeParameters
    MODULE PROCEDURE InitializeParameters
END INTERFACE


INTERFACE SimulationSelector
    MODULE PROCEDURE SimulationSelector
END INTERFACE 
!-------------------------------------------------------------------------------!
PUBLIC :: InitializeParameters
PUBLIC :: SimulationSelector
!-------------------------------------------------------------------------------!
! 
! 
! 
!===============================================================================!
CONTAINS
!===============================================================================!
! 
! 
! 
SUBROUTINE SimulationSelector()

USE MOD_FiniteVolume
USE MOD_MomentModels
USE MOD_TimeDiscretization
USE MOD_FiniteVolume_vars,ONLY: nGPsX, nGPsY, nDims, nXs, nYs, nMoms, nGhostsX, nGhostsY
USE MOD_FiniteVolume_vars,ONLY: ReconstrFlag, NCFluxFlag, ModelFlag, OutputFlag, PathFlag, TimeSchemeFlag, SplittingFlag

USE MOD_FiniteVolume_vars,ONLY: VarModlOpt, VarNConOpt, VarRecoOpt, VarOutsOpt, VarFluxOpt, VarOutsOpt, VarTimeOpt, VarTimeOpt
USE MOD_FiniteVolume_vars,ONLY: QuadPoints, QuadWeights, nGPath, QuadPathFlag, ViscosityFlag, RootsAH
USE MOD_FiniteVolume_vars,ONLY: QuadSourceWeights, QuadSourceNodes, nGSource
USE MOD_FiniteVolume_vars,ONLY: LegenPoly, DxLegenPoly
USE MOD_PhysicsFrame,ONLY: RootsLegendre01, LegendrePolyNormalized01
USE MOD_PhysicsFrame,ONLY: JacobiRoots

CHARACTER(LEN=255) :: ErrorMessage
REAL    :: xvals(1:nMoms), wwvals(1:nMoms)
INTEGER :: ii,jj
REAL    :: PPP, DPPP

PathFlag     = 1  ! Type of path, linear, quadratic or polynomial
QuadPathFlag = 3  ! Number of quadrature points for the path integral

IF (nMoms > 0) THEN
    CALL JacobiRoots(RootsAH)
END IF

SELECT CASE (ReconstrFlag)
CASE (1) ! NONE
    VarRecoOpt = "None, O(1)"
    nGhostsX = 1
    nGPsX    = 1 
CASE (2) ! MUSCL
    VarRecoOpt = "MUSCL, O(2)"
    nGhostsX = 1
    nGPsX    = 1
CASE (3) ! WENO3
    VarRecoOpt = "Weno3"
    nGhostsX = 1
    nGPsX    = 2
CASE (4) ! WENO5
    VarRecoOpt = "Weno5"
    nGhostsX = 2
    nGPsX    = 4
CASE DEFAULT
    ErrorMessage = "Reconstruction not implemented"
    WRITE (*, *) ErrorMessage
    STOP
END SELECT
#ifdef TWODIM
    nGhostsY = nGhostsX
    nGPsY    = nGPsX
#else
    nGhostsY = -1
    nGPsY    = 1
    nGPsX    = 1
#endif

!--------------------------------------------------------------------------------------!
! Path used to compute the integral corresponding to the linearzed Roe matrix
SELECT CASE (PathFlag)
CASE(1) ! Linear
    PathIntegral => LinearPath
CASE(2) ! Quadratic
    PathIntegral => QuadraticPath        
CASE(3) ! Power law negative
    PathIntegral => PowerlawMinusPath        
CASE(4) ! Power law positive
    PathIntegral => PowerlawPlusPath        
CASE DEFAULT
    PRINT *, "Path not implemented"
    STOP
END SELECT
!--------------------------------------------------------------------------------------!

SELECT CASE (ModelFlag)
CASE(1) ! "SWE1D"
    VarModlOpt = "SWE1D"
    SystemMatrix     => ZeroMatrix
    RemainMatrix     => ZeroMatrix
    ConservativeFlux => RiemannSolver
    SourceTerm       => SourceFuncSWME1D
CASE(2) ! "SWE2D"
    VarModlOpt = "SWE2D"
    SystemMatrix     => ZeroMatrix
    RemainMatrix     => ZeroMatrix
    ConservativeFlux => RiemannSolver
    SourceTerm       => ZeroSource    ! To set in an advance stage
CASE(3) ! "SWME1D"
    VarModlOpt = "SWME1D"
    SystemMatrix     => JacobianSWME1D
    RemainMatrix     => QmatrixSWME1D
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceFuncSWME1D
CASE(4) ! "HSWME1D"
    VarModlOpt = "HSWME1D"
    SystemMatrix     => JacobianHSWME1D
    RemainMatrix     => ZeroMatrix
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceFuncSWME1D
CASE(5) ! "BHSWME1D"
    VarModlOpt = "BetaHSWME1D"
    SystemMatrix     => JacobianBHSWME1D
    RemainMatrix     => ZeroMatrix
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceFuncSWME1D
CASE(6) ! "SWLME1D"
    VarModlOpt = "SWLME1D"
    SystemMatrix     => JacobianSWLME1D    
    RemainMatrix     => ZeroMatrix    
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceFuncSWME1D    
CASE(7) ! "ReducedModelI"
    VarModlOpt = "ReducedModelI"
    SystemMatrix     => JacobianReducedModelI  
    RemainMatrix     => ZeroMatrix    
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceReducedModelI

CASE(8) ! "HSWME1D-Inclined-NewtonianSlip"
    VarModlOpt = "HSWME1D-Inclined-NewtonianSlip"
    SystemMatrix     => JacobianInclinedHSWME1D
    RemainMatrix     => ZeroMatrix
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceNewtonianSlip

CASE(9) ! "Moment regularization models"
    VarModlOpt = "MHSWME1D-Moment"
    SystemMatrix     => JacobianMHSWME1D
    RemainMatrix     => ZeroMatrix
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceFuncSWME1D

CASE(10) ! "Moment regularization models"
    VarModlOpt = "PHSWME1D-Primitive"
    SystemMatrix     => JacobianPHSWME1D
    RemainMatrix     => ZeroMatrix
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceFuncSWME1D
    
CASE(11) ! "Moment regularization models"
    VarModlOpt = "PMHSWME1D-Moment-Primitive"
    SystemMatrix     => JacobianPHSWME1D
    RemainMatrix     => ZeroMatrix
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceFuncSWME1D
    
    !--------------------------------------------------------------------!   
CASE(12) ! "Source Newtonian"
    VarModlOpt = "HSWME-Newtonian-Slip"
    SystemMatrix     => JacobianInclinedHSWME1D
    RemainMatrix     => ZeroMatrix
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceNewtonianSlip     

CASE(13) ! "Source Newtonian Manning"
    VarModlOpt = "HSWME-Newtonian-Manning"
    SystemMatrix     => JacobianInclinedHSWME1D
    RemainMatrix     => ZeroMatrix
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceNewtonianManning     

CASE(14) ! "Source Savage Hutter"
    VarModlOpt = "HSWME-Savage-Hutter"
    SystemMatrix     => JacobianInclinedHSWME1D
    RemainMatrix     => ZeroMatrix
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceSavageHutter     

CASE(15) ! "Source Coulomb A"
    VarModlOpt = "HSWME-Coulomb-type-A"
    SystemMatrix     => JacobianInclinedHSWME1D
    RemainMatrix     => ZeroMatrix
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceCoulombA     
   
CASE(16) ! "Source Granular"
    VarModlOpt = "HSWME-Granular"
    SystemMatrix     => JacobianInclinedHSWME1D
    RemainMatrix     => ZeroMatrix
    ConservativeFlux => ZeroSolver
    SourceTerm       => SourceGranular           
    !--------------------------------------------------------------------!   
       
CASE(99) ! "SWME2D" ! Working on it
    VarModlOpt = "SWME2D"
    SystemMatrix     => JacobianSWME1D
    RemainMatrix     => QmatrixSWME1D
    ConservativeFlux => ZeroSolver
    SourceTerm       => ZeroSource ! To set in an advance stage
CASE DEFAULT
    PRINT *, "Error: 1D Model not implemented"
    STOP
END SELECT

PRINT*, achar(27)//'[1;96m    Running: "'//TRIM(VarModlOpt)//'" model'//achar(27)//'[0m'
PRINT*, achar(27)//'[1;96m'//REPEAT('-',50)//achar(27)//'[0m'

! Select time step:
SELECT CASE (ModelFlag)
CASE(1,2,99)
    TimeStep => TimeStepSWE
CASE(3)
    TimeStep => TimeStepGeneralSWME
CASE(4,5,8)
    TimeStep => TimeStepHSWME   
CASE(6)
    TimeStep => TimeStepSWLME
CASE DEFAULT
    TimeStep => TimeStepGeneralSWME
END SELECT

!--------------------------------------------------------------------------------------!
IF (SplittingFlag .EQ. 0) THEN
   SourceSplit => ZeroSource
ELSE
   SourceSplit => SourceTerm
   SourceTerm  => ZeroSource
END IF

!--------------------------------------------------------------------------------------!
SELECT CASE(NCFluxFlag)
CASE(1) ! PVM method
    NonConsMethod => PVM    
    SELECT CASE (ViscosityFlag)
    CASE(1) ! Lax-Friedrichs:
        VarNConOpt = "PVM-Lax-Friedrichs"
        ViscousMatrix => ViscousMatrixLxF
    CASE(2) ! Lax-Wendroff:
        VarNConOpt = "PVM-Lax-Wendroff"
        ViscousMatrix => ViscousMatrixLxW
    CASE(3) ! Force
        VarNConOpt = "PVM-Force"
        ViscousMatrix => ViscousMatrixForce
    CASE(4) ! HLL
        VarNConOpt = "PVM-HLL"
        ViscousMatrix => ViscousMatrixHLL     
    CASE DEFAULT
        PRINT *, achar(27)//'[93m    Error: PVM: Viscosity matrix not implemented'//achar(27)//'[0m'
        STOP
    END SELECT
CASE(2) ! Eigen Value decomposition method
    VarNConOpt = "Eigen-value-decomposition"
    NonConsMethod => EigenValueMethod
CASE DEFAULT
    PRINT *, achar(27)//'[93m    Error: Non-conservative method not implemented'//achar(27)//'[0m'
    STOP
END SELECT
!--------------------------------------------------------------------------------------!

SELECT CASE (OutputFlag)
CASE(1)
    VarOutsOpt = "OCTAVE"
CASE(2)
    VarOutsOpt = "TECPLOT"
CASE(3)
    VarOutsOpt = "OCTAVE & TECPLOT"
CASE(4)
    VarOutsOpt = "New"
END SELECT

VarFluxOpt = "Not in use"
!VarFluxOpt = "Russanov"


!-------------------------------------------------------------------------------!
SELECT CASE (TimeSchemeFlag)
CASE (1)
    VarTimeOpt = "Forward-Euler"
    TimeApproximation => TimeDiscretizationByForwardEuler
CASE (2)
    VarTimeOpt = "SSPRK4"
    TimeApproximation => TimeDiscretizationBySSPRK4
CASE (3)
    VarTimeOpt = "RK65"
    TimeApproximation => TimeDiscretizationByRK65
CASE (4)
    VarTimeOpt = "DeC5"
    TimeApproximation => TimeDiscretizationByDeC5
CASE DEFAULT
    PRINT *, '\n',achar(27)//'[1;103;31m'//REPEAT(' ',4)//'Error: ¡Time discretization not implemented!   '//achar(27)//'[0m'
    STOP
END SELECT
!-------------------------------------------------------------------------------!
! Gauss quadrature - Legendre polynomials in [-1,1]:
SELECT CASE(QuadPathFlag) 
CASE(1) 
    nGPath = 1
    QuadPoints(1)  = 0.0
    QuadWeights(1) = 2.0
CASE(2)
    nGPath = 2
    QuadPoints(1:2)  = (1.0/sqrt(3.))*(/-1.0, 1.0/)
    QuadWeights(1:2) = (/1.0, 1.0/)
CASE(3)
    nGPath = 3
    QuadWeights(1:3) = (  1.0/9.0)*(/  5.0,  8.0,  5.0 /)
    QuadPoints(1:3)  = sqrt(3./5.)*(/ -1.0,  0.0,  1.0 /)
CASE(4)
    nGPath = 4
    QuadWeights(1:4) = (1./36.)*(18.0*(/1.,1.,1.,1./)  +  SQRT(30.)*(/1.0,-1.0,-1.0,1.0 /))
    QuadPoints(1:4)  = (/-SQRT(3./7. + 2./7.*SQRT(6./5.)), -SQRT(3./7. - 2./7.*SQRT(6./5.)),  &
                          SQRT(3./7. - 2./7.*SQRT(6./5.)),  SQRT(3./7. + 2./7.*SQRT(6./5.))   /)                        

END SELECT
! Gauss Quadrature rules for integral between 0 and 1;  sum 0.5*w_i* f( x_i + 0.5 )
! Here the 0.5 everywhere is due to the (b-a)/2 factor in the quadratures w_i, x_i from the interval [-1,1]
QuadPoints  = 0.5*(QuadPoints + 1.0) 
QuadWeights = 0.5*QuadWeights

!---------------------------------------------------------------------------------------------------------------------
! Quadrature rule used only for the source terms:

nGSource = 8
ALLOCATE(QuadSourceWeights(1:nGSource), QuadSourceNodes(1:nGSource), & 
         LegenPoly(1:nGSource,1:nMoms), DxLegenPoly(1:nGSource,1:nMoms))

QuadSourceWeights = 0.0; 
QuadSourceNodes   = 0.0;
LegenPoly         = 0.0; 
DxLegenPoly       = 0.0;

SELECT CASE(nGSource)
CASE(4)    
    QuadSourceWeights = (1./36.)*( 18.*(/1.,1.,1.,1./) & 
                           + SQRT(30.)*(/1.,-1.,-1.,1./)  )
                        
    QuadSourceNodes   = 0.5*(/ -SQRT(3./7. + 2./7.*SQRT(6./5.)), &
                               -SQRT(3./7. - 2./7.*SQRT(6./5.)), &
                                SQRT(3./7. - 2./7.*SQRT(6./5.)), &
                                SQRT(3./7. + 2./7.*SQRT(6./5.))/)

CASE(5)
    QuadSourceWeights = (1./900.)*(322.*(/1.,1., (128./255.)*(900./322.),1.,1./) &
                        + 13*SQRT(70.0)*(/-1.,1.,0.,1.,-1./) )
     
    QuadSourceNodes   =  (1.0/3.0)*(/ -SQRT(5.0 + 2.0*SQRT(10.0/7.0)), &
                                      -SQRT(5.0 - 2.0*SQRT(10.0/7.0)), &
                                       0.0, &
                                       SQRT(5.0 - 2.0*SQRT(10.0/7.0)), &                                         
                                       SQRT(5.0 + 2.0*SQRT(10.0/7.0)) /)  
                                       
CASE(6)

    QuadSourceWeights = (/0.17132449237917034500, &
                          0.36076157304813860757, &
                          0.46791393457269104740, &
                          0.46791393457269104740, & 
                          0.36076157304813860757, & 
                          0.17132449237917034500 /)
    
    QuadSourceNodes   = (/-0.93246951420315202780, &
                          -0.66120938646626451366, &
                          -0.23861918608319690863, &
                           0.23861918608319690863, &
                           0.66120938646626451366, &
                           0.93246951420315202780/)

CASE(7)
    QuadSourceWeights = (/0.12948496616886969327, & 
                          0.27970539148927666790, &
                          0.38183005050511894495, &
                          0.41795918367346938776, &
                          0.38183005050511894495, &
                          0.27970539148927666790, &
                          0.12948496616886969327 /)                          

    QuadSourceNodes   = (/-0.94910791234275852453, &
                          -0.74153118559939443986, &
                          -0.40584515137739716691, &
                           0.00000000000000000000, &
                           0.40584515137739716691, &
                           0.74153118559939443986, &
                           0.94910791234275852453/)  
CASE(8)
    QuadSourceWeights = (/-0.96028985649753623168, &
                          -0.79666647741362673959, &
                          -0.52553240991632898582, &
                          -0.18343464249564980494, &
                           0.18343464249564980494, &
                           0.52553240991632898582, &
                           0.79666647741362673959, &
                           0.96028985649753623168/)

    QuadSourceWeights = (/ 0.10122853629037625915, &
                           0.22238103445337447054, &
                           0.31370664587788728734, &
                           0.36268378337836198297, &
                           0.36268378337836198297, &
                           0.31370664587788728734, &
                           0.22238103445337447054, &
                           0.10122853629037625915/)                                   
END SELECT
! convert to the interval [0,1]
QuadSourceWeights = 0.5*QuadSourceWeights
QuadSourceNodes   = 0.5*(QuadSourceNodes + 1.0)
!---------------------------------------------------------------------------------------------------------------------
! General quadrature rule with (w_i, xi_i) for Legendre polynomial [-1,1]
!
! int_a^b f(x) dx = ((b-a)/2) * sum w_i * f( (b-a)/2 *xi_i + (a_i+b_i)/2) 
!
! when a = 0 and b = 1 then
!      int_0^1 f(x) dx  =  sum  0.5*w_i * f(  0.5*(xi_i + 1.0)  ) 
DO ii = 1, nGSource 
   DO jj = 1, nMoms   
       ! These are Phi_j(x_i)  in the order (i,j)   
       CALL LegendrePolyNormalized01(QuadSourceWeights(ii), LegenPoly(ii,jj), DxLegenPoly(ii,jj), jj)
   END DO
END DO

END SUBROUTINE SimulationSelector




!===============================================================================!
SUBROUTINE InitializeParameters()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: PI, CFL, TEnd, Gravity
USE MOD_FiniteVolume_vars,ONLY: nXs, nYs, MESH_X0, MESH_X1
USE MOD_FiniteVolume_vars,ONLY: ReconstrFlag, ReconstrFixFlag, BoundaryFlag
USE MOD_FiniteVolume_vars,ONLY: OutputFlag, InitialFlag, nOutputFiles, TimeSchemeFlag
USE MOD_FiniteVolume_vars,ONLY: VarNameVisu, BathymetryFlag, ViscosityFlag
!-------------------------------------------------------------------------------!
! Moment related matrices
USE MOD_FiniteVolume_vars,ONLY: Acoef, Bcoef, Ccoef
USE MOD_FiniteVolume_vars,ONLY: SubIdentity, Identity
USE MOD_FiniteVolume_vars,ONLY: VecDj, VectorBeta
USE MOD_FiniteVolume_vars,ONLY: ParameterNu, ParameterLambda
USE MOD_FiniteVolume_vars,ONLY: nMoms, nVar
USE MOD_FiniteVolume_vars,ONLY: ModelFlag, NCFluxFlag, SplittingFlag, Examplefolder
USE MOD_FiniteVolume_vars,ONLY: QuadPoints, QuadWeights, AngleIncl, AuxiliaryParameter



!------------
USE MOD_FiniteVolume_vars,ONLY: nMomsOut, RootsAH

IMPLICIT NONE
!-------------------------------------------------------------------------------!
CHARACTER(LEN=255) :: ErrorMessage, auxarg, Auxfoldername
INTEGER            :: iarg, nargs  
CHARACTER(len=255) :: arg
INTEGER            :: ii, jj, kk, ndirs, nauxdirs, markdir

REAL :: newc(11,11)

newc = 0.0
!-------------------------------------------------------------------------------!
nauxdirs = 0; ndirs = 0; markdir  = 0 

QuadPoints  = 0.0 
QuadWeights = 0.0

!====================================================================================!
! PREAMBLE SECTION TO CHECK FOLDER AND INPUT FILE, AND READ FROM INPUT
!------------------------------------------------------------------------------------!
nargs = command_argument_COUNT()
!call get_command_argument_COUNT(nargs)
print*, "number of args:", nargs
IF (nargs > 0) THEN 
    CALL get_command_ARGUMENT(1, arg)
    Examplefolder = trim(arg)
END IF

CALL execute_command_line("find "//Examplefolder//" -type d > dummy.txt" ) ! dummy file to confirm folder
OPEN(unit = 99, file = 'dummy.txt')
READ(99,*, iostat = nauxdirs)
IF (nauxdirs/=0) THEN; 
    PRINT*,"====================\nThe folder is not found, the program stops here.\n===================="
    STOP
END IF
nauxdirs = 0
OPEN(unit = 88, file = Examplefolder//"/input.txt")
READ(88,*, iostat = nauxdirs)
IF (nauxdirs/=0) THEN; 
    PRINT*, "====================\nThe program requires an input.txt file to be executed\n===================="
    STOP
END IF
CLOSE(88)
CLOSE(99)
CALL execute_command_line("mkdir -p "//Examplefolder//"/Error")  ! Creates Error folder
CALL execute_command_line("rm "//Examplefolder//"/*.dat" )       ! Clean the dat files
CALL execute_command_line("rm "//Examplefolder//"/**/*.dat" )    ! Clean the dat files
OPEN(UNIT = 1, FILE = Examplefolder//'/input.txt', STATUS = 'OLD', ACTION = 'READ')
READ(1,*) ModelFlag
READ(1,*) nVar
READ(1,*) nMoms
READ(1,*) nMomsOut
READ(1,*) TEnd
READ(1,*) Gravity
READ(1,*) ParameterNu
READ(1,*) ParameterLambda
READ(1,*) nXs
READ(1,*) nYs
READ(1,*) MESH_X0(1) ! left
READ(1,*) MESH_X1(1) ! right
READ(1,*) MESH_X0(2) ! bottom 
READ(1,*) MESH_X1(2) ! top
READ(1,*) BoundaryFlag(1) ! left
READ(1,*) BoundaryFlag(2) ! right
READ(1,*) BoundaryFlag(3) ! bottom
READ(1,*) BoundaryFlag(4) ! top
READ(1,*) BathymetryFlag
!----------------------------
READ(1,*) InitialFlag
READ(1,*) ReconstrFlag 
READ(1,*) TimeSchemeFlag
READ(1,*) OutputFlag
READ(1,*) nOutputFiles
READ(1,*) CFL
READ(1,*) NCFluxFlag
READ(1,*) ViscosityFlag
READ(1,*) SplittingFlag
READ(1,*) AngleIncl
READ(1,*) AuxiliaryParameter
!----------------------------
CLOSE (1)

! The AuxiliaryParameter is used in case of needed in some models, in particular in the source terms
!nVar  = 1 + nDims + nMoms*nDims
AngleIncl = AngleIncl*PI

!====================================================================================!
! MOMENT MODEL MATRICES (LEGENDRE POLYNOMIAL DEPENDENT)
!------------------------------------------------------------------------------------!

ALLOCATE(Identity(1:nvar,1:nvar), SubIdentity(1:max(1,nMoms), 1:max(1,nMoms)), &
         VecDj(1:max(1,nMoms)),   VectorBeta(1:nVar), VarNameVisu(1:nVar+1), Ccoef(1:nMoms,1:nMoms), RootsAH(1:nMoms))

SubIdentity = 0.0 ! N x N
Identity    = 0.0 ! nVar x nVar
VecDj       = 0.0 
VectorBeta  = 0.0
Acoef       = 0.0
Bcoef       = 0.0
Ccoef       = 0.0

DO jj = 1, nMoms
   VecDj(jj)           = 1./(2.*REAL(jj) + 1.0)
   SubIdentity(jj, jj) = 1.0
END DO
DO jj = 1, nVar
   VectorBeta(jj)   = 0.0
   Identity(jj, jj) = 1.0
END DO
!-------------------------------------------------------------------------------------!
Acoef(1,:,:) = transpose(reshape((/               &
                0.,  2./5., 0., 0., 0.,           &
                2./5., 0., 9./35., 0., 0.,        &  
                0., 9./35., 0., 4./21., 0.,       &
                0., 0., 4./21., 0., 5./33.,       &
                0., 0., 0., 5./33., 0./),         &
                (/5, 5/)))
Acoef(2,:,:) = transpose(reshape((/               &
                2./3., 0., 3./7., 0., 0.,         &
                0., 2./7., 0., 2./7., 0.,         &
                3./7., 0., 4./21., 0., 50./231.,  &
                0., 2./7., 0., 100./693., 0.,     &
                0., 0., 50./231., 0., 50./429./), &
                (/5, 5/)))
Acoef(3,:,:) = transpose(reshape((/               &
                0., 3./5., 0., 4./9., 0.,         &
                3./5., 0., 4./15., 0., 10./33.,   &
                0., 4./15., 0., 2./11., 0.,       &
                4./9., 0., 2./11., 0., 20./143.,  &
                0., 10./33., 0., 20./143., 0./),  &
                (/5, 5/)))
Acoef(4,:,:) = transpose(reshape((/                 &
                0., 0., 4./7., 0., 5./11.,          &
                0., 18./35., 0., 20./77., 0.,       &
                4./7., 0., 18./77., 0., 180./1001., &
                0., 20./77., 0., 162./1001., 0.,    &
                5./11., 0., 180./1001., 0., 18./143./), &
                (/5, 5/)))
Acoef(5,:,:) = transpose(reshape((/               &   
                0., 0., 0., 5./9., 0.,            &
                0., 0., 10./21., 0., 10./39.,     &
                0., 10./21., 0., 20./91., 0.,     &
                5./9., 0., 20./91., 0., 2./13.,   &
                0., 10./39., 0., 2./13., 0. /),   &
                (/5, 5/)))
!==================================================
Bcoef(1,:,:) = transpose(reshape((/               &   
                0., 1./5., 0., 0., 0.,            &
               -1./5., 0., 3./35., 0., 0.,        &
                0., -3./35., 0., 1./21., 0.,      &
                0., 0., -1./21., 0., 1./33.,      &
                0., 0., 0., -1./33., 0./),        &
                (/5, 5/)))
Bcoef(2,:,:) = transpose(reshape((/               &
               -1., 0., 3./7., 0., 0.,            &
                0., -1./7., 0., 4./21., 0.,       &
               -2./7., 0., -1./21., 0., 25./231., &
                0., -1./7., 0., -5./231., 0.,     &
                0., 0., -20./231., 0.,-5./429./), &
                (/5, 5/)))
Bcoef(3,:,:) = transpose(reshape((/               &
                0., -6./5., 0., 2./3., 0.,        &
               -4./5., 0., -2./15., 0., 10./33.,  &
                0., -1./5., 0., -1./33., 0.,      &
               -1./3., 0., -1./11., 0., -1./143., &
                0., -2./11., 0., -2./39., 0. /),  &
                (/5, 5/)))
Bcoef(4,:,:) = transpose(reshape((/               &
                0., 0., -10./7., 0., 10./11.,     &
                0., -6./7., 0., -10./77., 0.,     &
               -5./7., 0., -15./77., 0., -15./1001.,     &
                0., -17./77., 0., -81./1001., 0.,        &
               -4./11., 0., -114./1001., 0., -6./143./), &
                (/5, 5/)))
Bcoef(5,:,:) = transpose(reshape((/               &
                0., 0., 0., -5./3., 0.,           &
                0., 0., -20./21., 0., -5./39.,    &
                0., -5./7., 0., -55./273., 0.,    &
               -2./3., 0., -19./91., 0., -1./13., &
                0., -3./13., 0., -4./39., 0.  /), &
                (/5, 5/)))
!==================================================
! Ccoef(:,:) = transpose(reshape((/      &
!               4.,  0.,  4.,  0.,  4.,      &
!               0., 12.,  0., 12.,  0.,    &
!               4.,  0., 24.,  0., 24.,    &
!               0., 12.,  0., 40.,  0.,    &
!               4.,  0., 24.,  0., 60./),  &
!               (/5, 5/)))
               
DO ii = 2, nMoms+1
    DO jj = 1, nMoms
        IF (mod(ii + jj, 2) == 0) THEN
            Ccoef(ii-1,jj) = 0.0
        ELSE
            Ccoef(ii-1,jj) = 0.5 * REAL(MIN(ii - 1, jj)) * REAL(MIN(ii - 1, jj) + 1)
        END IF
    END DO
END DO               
Ccoef = 4.0*Ccoef

!====================================================================================!
! TITLE CONTAINING THE UNKNOWN VARIABLES
VarNameVisu(1)      = "Depth"
VarNameVisu(2)      = "VelocityX"
VarNameVisu(nVar+1) = "Bathymetry"
VarNameVisu(3)      = "VelocityY" ! Overwritten when 1D
#ifndef TWODIM
nYs = 1  ! Modify the matrix size corresponding to the y-axis when 1D
DO ii = 1, nMoms
    WRITE(VarNameVisu(ii + 2),'(I0)') ii
    VarNameVisu(ii + 2) = "MomentX_i"//VarNameVisu(ii)
END DO
#else
DO ii = 1, nMoms
    WRITE(VarNameVisu(2*ii-1 + 3),'(I0)') ii
    WRITE(VarNameVisu(2*ii   + 3),'(I0)') ii
    VarNameVisu(2*ii-1 + 3) = "MomentX_"//VarNameVisu(2*ii-1 + 3)
    VarNameVisu(2*ii   + 3) = "MomentY_"//VarNameVisu(2*ii   + 3)
END DO
#endif
!================================================================================ VarInitOpt = "Wave-dry-island" ====!
! SHORT PRINTING SECTION
PRINT*, "\n ", REPEAT('=', 50)
#ifdef TWODIM
    PRINT*, "2D simulation"
#else
    PRINT *, achar(27)//'[1;104m'//'1D simulation'//repeat(' ',37)//achar(27)//'[0m'
#endif
PRINT*, REPEAT('=', 50)
#ifdef WELLBALANCED
    PRINT*, "Well balanced"
#endif
!-------------------------------------------------------------------------------!
PRINT *, achar(27)//'[1;106m'//'Basic data'//repeat(' ',40)//achar(27)//'[0m'
PRINT*, REPEAT('-', 50)
PRINT('(A,I0)'), " nXs               = ", nXs
PRINT('(A,I0)'), " nYs               = ", nYs
PRINT('(A,I0)'), " Test -Initial C.  = ", InitialFlag
PRINT('(A,I0)'), " Time Scheme       = ", TimeSchemeFlag
PRINT('(A,I0)'), " Reconstruction    = ", ReconstrFlag
PRINT('(A,I0)'), " ReconstructionFix = ", ReconstrFixFlag
PRINT*, REPEAT('=', 50)
PRINT('(A)'), " This program currently supports N < 6 moments"
PRINT('(A)'), " (Add bigger matrices to parameters.f90)"
PRINT*, achar(27)//'[1;97m'//REPEAT('-',50)//achar(27)//'[0m'
ReconstrFixFlag = ReconstrFlag

!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeParameters
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
