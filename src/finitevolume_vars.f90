!===============================================================================!
MODULE MOD_FiniteVolume_vars
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PUBLIC
!-------------------------------------------------------------------------------!
! >> GLOBAL VARIABLES                                                           !
!-------------------------------------------------------------------------------!
#ifdef TWODIM
INTEGER,PARAMETER   :: nDims = 2 ! Two dimensional
#else
INTEGER,PARAMETER   :: nDims = 1 ! One dimensional
#endif
INTEGER             :: nVar          ! Total number of unknowns
INTEGER             :: nMoms         ! Number of Moments within vector solution
INTEGER             :: nMomsOut      ! Additional moments, used only in a post-processing routine

REAL                :: MESH_SX(1:2)  ! 
REAL                :: MESH_X0(1:2)  ! Pair (a,c), if the domain is [a,b]x[c,d]
REAL                :: MESH_X1(1:2)  ! Pair (b,d), if the domain is [a,b]x[c,d]
REAL,ALLOCATABLE    :: MESH_DX(:)    ! Mesh sizes (dx,dy) = ((b-a)/nXs, (d-c)/nYs)
INTEGER             :: nGPsX         ! Gauss-Quadrature points, reconstruction in x-direction
INTEGER             :: nGPsY         ! Gauss-Quadrature points, reconstruction in y-direction
INTEGER             :: nGhostsX      ! Number of Ghost cells in x-direction
INTEGER             :: nGhostsY      ! Number of Ghost cells in y-direction
INTEGER             :: nXs           ! Number of cells in x-direction (within the domain)
INTEGER             :: nYs           ! Number of cells in y-direction (within the domain)

REAL,ALLOCATABLE    :: MeshNodes(:,:,:)    ! Set the nodes of the mesh
REAL,ALLOCATABLE    :: MeshBary(:,:,:)
REAL,ALLOCATABLE    :: MeshGP(:,:,:,:,:)
REAL,ALLOCATABLE    :: WeightsGP(:,:)
REAL,ALLOCATABLE    :: WeightsGPBnd(:)
REAL,ALLOCATABLE    :: NormVectX(:,:,:,:)  ! Normal vector at the edge, x-component
REAL,ALLOCATABLE    :: NormVectY(:,:,:,:)  ! Normal vector at the edge, y-component
REAL,ALLOCATABLE    :: TangVectX(:,:,:,:)  ! Tangent vector at the edge, x-component
REAL,ALLOCATABLE    :: TangVectY(:,:,:,:)  ! Tangent vector at the edge, y-component

!-------------------------------------------------------------------------------------!
! Unknowns:
REAL,ALLOCATABLE    :: U(:,:,:)       ! Vector of unknowns: Conservative
REAL,ALLOCATABLE    :: V(:,:,:)       ! Vector of unknowns: Primitive
!-------------------------------------------------------------------------------------!
REAL,ALLOCATABLE    :: SSFriction(:,:,:)  ! Source vector
REAL,ALLOCATABLE    :: SWB(:,:,:)
REAL,ALLOCATABLE    :: Ut(:,:,:)
REAL,ALLOCATABLE    :: FX(:,:,:)      ! Total numerical flux in the x-direction
REAL,ALLOCATABLE    :: FY(:,:,:)      ! Total numerical flux in the y-direction
REAL,ALLOCATABLE    :: FXWB(:,:,:)
REAL,ALLOCATABLE    :: FYWB(:,:,:)
REAL,ALLOCATABLE    :: WM(:,:,:,:)    ! Polynomial reconstruction "minus"
REAL,ALLOCATABLE    :: WP(:,:,:,:)    ! Polynomial reconstruction "plus"
REAL,ALLOCATABLE    :: AuxFlux(:,:)      ! Auxiliary numerical flux at quadrature points
LOGICAL,ALLOCATABLE :: Ind(:,:,:)     ! To save the indexes from the shock indicator

REAL,ALLOCATABLE    :: UN0(:,:,:)     ! Auxiliary vectors, time discretization, RK methods
REAL,ALLOCATABLE    :: K0(:,:,:)      ! """
REAL,ALLOCATABLE    :: K1(:,:,:)      ! """ 
REAL,ALLOCATABLE    :: K2(:,:,:)      ! """
REAL,ALLOCATABLE    :: K3(:,:,:)      ! """
REAL,ALLOCATABLE    :: K4(:,:,:)      ! """
REAL,ALLOCATABLE    :: K5(:,:,:)      ! """

REAL,ALLOCATABLE    :: FUp(:,:,:,:)   ! Auxiliary vectors, used only in the DeC5 method
REAL,ALLOCATABLE    :: Up(:,:,:,:)    ! """
REAL,ALLOCATABLE    :: Ua(:,:,:,:)    ! """

REAL,ALLOCATABLE    :: UtWB(:,:,:)
REAL,ALLOCATABLE    :: Bath(:,:)      ! Cell average of the bathymetry at quadrature points
REAL,ALLOCATABLE    :: dBath(:,:)     ! Cell average of the derivative of the bathymetry at quadrature points
REAL,ALLOCATABLE    :: BvecX(:)       ! Pointwise evaluation of the bathymetry at cell centers plus borders
REAL,ALLOCATABLE    :: BvecXrot(:)    ! Auxiliary

REAL,ALLOCATABLE    :: AuxDplus(:)    ! Fluxtuation plus, auxiliary variable
REAL,ALLOCATABLE    :: AuxDmins(:)    ! Fluxtuation plus, auxiliary variable
REAL,ALLOCATABLE    :: DplusX(:,:,:)  ! Total fluxtuation plus, x-direction
REAL,ALLOCATABLE    :: DplusY(:,:,:)  ! Total fluxtuation plus, y-direction
REAL,ALLOCATABLE    :: DminsX(:,:,:)  ! Total fluxtuation minus, x-direction
REAL,ALLOCATABLE    :: DminsY(:,:,:)  ! Total fluxtuation minus, y-direction


REAL,ALLOCATABLE    :: PrimRefState1(:)
REAL,ALLOCATABLE    :: PrimRefState2(:)
REAL,ALLOCATABLE    :: PrimRefState3(:)
REAL,ALLOCATABLE    :: PrimRefState4(:)

INTEGER,PARAMETER   :: UNIT_FILE = 123
INTEGER             :: nOutputFiles

REAL                :: t          ! time
REAL                :: tGlobal    ! 
REAL                :: dt         ! delta t, time step
REAL                :: dt_save    ! time step to save
REAL                :: CFL        ! constant used to reduce the CFL condition
REAL                :: tEnd       ! end time point
REAL                :: Gravity    ! constant of gravity
REAL                :: LambdaMaxX ! maximum eigenvalue, x-direction
REAL                :: LambdaMaxY ! maximum eigenvalue, y-direction

REAL,PARAMETER      :: wLobatto = 1./12.   ! weight at Gauss-Lobatto scheme
REAL,PARAMETER      :: WENOEPS  = 1.0E-24  ! Epsilon at WENO scheme
INTEGER,PARAMETER   :: WENOEXP  = 2.0      ! Power at WENO scheme

REAL,PARAMETER      :: PI           = ACOS(-1.0)
REAL,PARAMETER      :: EPS          = 1.0E-6
REAL,PARAMETER      :: ACCURACY     = 1.0E-30
REAL,PARAMETER      :: MIN_SPEED    = 0.0
REAL,PARAMETER      :: MIN_TIMESTEP = 1.0E-30
REAL,PARAMETER      :: Hmin = 1.0E-6

CHARACTER(LEN=255),ALLOCATABLE :: VarNameVisu(:)

!-------------------------------------------------------------------------------------!
! Moment model equations:nMoms
REAL    :: ParameterNu
REAL    :: ParameterLambda
REAL    :: Acoef(1:5, 1:5, 1:5)
REAL    :: Bcoef(1:5, 1:5, 1:5)

REAL, ALLOCATABLE   :: Ccoef(:,:)       ! Matrix C_ij
REAL, ALLOCATABLE   :: Identity(:,:)    ! Identity matrix of size (nVar+2)*(nVar+2), 1D models only
REAL, ALLOCATABLE   :: SubIdentity(:,:) ! Identity matrix of size nVar*nVar, 1D models only
REAL, ALLOCATABLE   :: VecDj(:)         ! Auxiliary vector 1/(2j+1)
REAL, ALLOCATABLE   :: VectorBeta(:)    ! Beta vector used at the Beta-HSWME

REAL, ALLOCATABLE   :: RootsAH(:)       ! Roots characteristic polynomial arising in the A_H system
!-------------------------------------------------------------------------------------!
!Flags
INTEGER :: BathymetryFlag    ! Selects the bathymetry
INTEGER :: ReconstrFlag      ! Selects the polynomial reconstruction
INTEGER :: ReconstrFixFlag   ! TODO: remove this part.
INTEGER :: TimeSchemeFlag    ! Selects the time approximation 
INTEGER :: SplittingFlag     ! Selects the splitting treatment of the source term
INTEGER :: ModelFlag         ! Selects the model (matrices, fluxes and source function)
INTEGER :: NCFluxFlag        ! Selects the non-conservative method
INTEGER :: ViscosityFlag     ! Selects the viscosity matrix from the PVM method
INTEGER :: PathFlag          ! Selects the path at the Generalized Roe matrix
INTEGER :: OutputFlag        ! Selects the style of the output file
INTEGER :: InitialFlag       ! Selects the initial condition
INTEGER :: QuadPathFlag      ! Selects the order of the quadrature rule for the path integrals
INTEGER :: BoundaryFlag(4)   ! Selects the type of each boundary condition
!-------------------------------------------------------------------------------------!
! Labels 
CHARACTER(LEN=255) :: VarTimeOpt    ! Time Approximation
CHARACTER(LEN=255) :: VarNConOpt    ! Non-conservative approximation
CHARACTER(LEN=255) :: VarFluxOpt    ! Numerical Flux
CHARACTER(LEN=255) :: VarRecoOpt    ! Polynomial reconstruction
CHARACTER(LEN=255) :: VarOutsOpt    ! Output format
CHARACTER(LEN=255) :: VarInitOpt    ! Initial Condition 
CHARACTER(LEN=255) :: VarBConOpt(4) ! Boundary conditions
CHARACTER(LEN=255) :: VarModlOpt    ! Models
CHARACTER(LEN=255) :: VarBathOpt    ! Bathymetry
CHARACTER(LEN=255) :: VarFricOpt    ! Friction 
CHARACTER(LEN=255) :: auxBC         
CHARACTER(len=:), ALLOCATABLE :: Examplefolder  ! Folder name of the current example

!-------------------------------------------------------------------------------------!
! Here we set the length 4, increase in case of more quadrature points.
REAL    :: QuadPoints(1:4)   ! Store the quadrature points for the path integrals
REAL    :: QuadWeights(1:4)  ! Store the weights for the path integrals
INTEGER :: nGPath            ! Number of quadrature points used (depends on PathFlag)

!-------------------------------------------------------------------------------------!
INTEGER :: IndexPass(2)         ! Auxiliar 2-vector to store indexes for the source terms
REAL, ALLOCATABLE   :: NNv(:,:) ! To pass the normal vector
REAL, ALLOCATABLE   :: TTv(:,:) ! To pass the tangent vector

!-------------------------------------------------------------------------------------!
! Savage-Hutter Moment Equations
REAL :: AngleIncl                 ! pi/6, angle used for the inclined model
REAL,PARAMETER :: EpsilonRatio = 0.01        ! This is H/L to properly define later on
REAL,PARAMETER :: ParameterViscosity = 0.1   ! Viscosity mu

!-----------------------------------
REAL,ALLOCATABLE     :: QuadSourceWeights(:) ! Fourth order quadrature rule used for the source terms
REAL,ALLOCATABLE     :: QuadSourceNodes(:)   ! Fourth order quadrature rule used for the source terms
INTEGER              :: nGSource             ! 
REAL,ALLOCATABLE     :: LegenPoly(:,:), DxLegenPoly(:,:)
!-----------------------------------

REAL :: AuxiliaryParameter

END MODULE MOD_FiniteVolume_vars
!===============================================================================!
