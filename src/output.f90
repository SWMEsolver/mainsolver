!===============================================================================!
MODULE MOD_Output
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE WriteMeshToDisk
    MODULE PROCEDURE WriteMeshToDisk
END INTERFACE

INTERFACE WriteSolutionToDisk
    MODULE PROCEDURE WriteSolutionToDisk
END INTERFACE

INTERFACE ReducedModelPostProcessing
    MODULE PROCEDURE ReducedModelPostProcessing
END INTERFACE

INTERFACE ErrorCalculation  
  MODULE PROCEDURE ErrorCalculation 
END INTERFACE

!-------------------------------------------------------------------------------!
PUBLIC :: WriteMeshToDisk
PUBLIC :: WriteSolutionToDisk
PUBLIC :: ReducedModelPostProcessing
PUBLIC :: ErrorCalculation
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
!===============================================================================!
SUBROUTINE WriteSolutionToDisk(tOut,iter)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nVar, nDims, nXs, nYs
USE MOD_FiniteVolume_vars,ONLY: nGPsX, nGPsY
USE MOD_FiniteVolume_vars,ONLY: U, V
USE MOD_FiniteVolume_vars,ONLY: t, tGlobal
USE MOD_FiniteVolume_vars,ONLY: MeshBary, MeshGP, Bath, WeightsGP
USE MOD_FiniteVolume_vars,ONLY: VarNameVisu, Examplefolder, OutputFlag
USE MOD_FiniteVolume_vars,ONLY: InitialFlag
USE MOD_PhysicsFrame,           ONLY: ExactFunction
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: tOut
INTEGER,INTENT(IN) :: iter
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj, iGP, jGP
REAL,ALLOCATABLE   :: U_NVisu(:,:,:)
REAL,ALLOCATABLE   :: Error(:,:,:)
REAL,ALLOCATABLE   :: Uexact(:,:,:)
CHARACTER(LEN=255) :: FileNameTecU, FileNameOctU, FileNameNewU
CHARACTER(LEN=255) :: FileNameTecE, FileNameOctE, FileNameNewE
CHARACTER(LEN=4)   :: FVDim 
REAL, DIMENSION(nVar, nGPsX, nGPsY) :: Utemp
CHARACTER(LEN=255) :: TimeStamp
!-------------------------------------------------------------------------------!
#ifdef TWODIM
    FVDim = "FV2D"
#else
    FVDim = "FV1D"
#endif
!-------------------------------------------------------------------------------!
ALLOCATE(U_NVisu(1:nVar + 1,  1:nXs,  1:nYs))
ALLOCATE(Error(  1:nVar + 1,  1:nXs,  1:nYs))
ALLOCATE(Uexact( 1:nVar,      1:nXs,  1:nYs))
!-------------------------------------------------------------------------------!
WRITE(TimeStamp,'(I0)') iter
FileNameOctU = "soln"//TRIM(TimeStamp)
FileNameTecU = "soln"//TRIM(TimeStamp)
FileNameNewU = "soln"//TRIM(TimeStamp)
FileNameOctE = "error"//TRIM(TimeStamp)
FileNameTecE = "error"//TRIM(TimeStamp)
FileNameNewE = "error"//TRIM(TimeStamp)
!-------------------------------------------------------------------------------!
FileNameOctU = TRIM(Examplefolder//"/"//TRIM(FileNameOctU)//"_"//FVDim(3:4))
FileNameTecU = TRIM(Examplefolder//"/"//TRIM(FileNameTecU)//"_"//FVDim(3:4))
FileNameNewU = TRIM(Examplefolder//"/"//TRIM(FileNameNewU)//"_"//FVDim(3:4))
FileNameOctE = TRIM(Examplefolder//"/Error/"//TRIM(FileNameOctE)//"_"//FVDim(3:4))
FileNameTecE = TRIM(Examplefolder//"/Error/"//TRIM(FileNameTecE)//"_"//FVDim(3:4))
FileNameNewE = TRIM(Examplefolder//"/Error/"//TRIM(FileNameNewE)//"_"//FVDim(3:4))
!-------------------------------------------------------------------------------!
! Error calculation with the exact solution at each cell:
U_NVisu = 0.
Uexact  = 0.
Utemp   = 0.

DO jj = 1, nYs
    DO ii = 1, nXs
        U_NVisu(1:nVar, ii, jj)   = V(1:nVar, ii, jj)
        U_NVisu(nVar + 1, ii, jj) = Bath(ii, jj)
        DO iGP = 1, nGPsX
            DO jGP = 1, nGPsY
                CALL ExactFunction(tGlobal, MeshGP(:, ii, jj, iGP, jGP), Utemp(:, iGP, jGP))
                Uexact(:, ii, jj) = Uexact(:, ii, jj) + WeightsGP(iGP, jGP)*Utemp(:, iGP, jGP)
            END DO
        END DO
        Error(1:nVar, ii, jj)   = U(1:nVar, ii, jj) - Uexact(1:nVar, ii, jj)
        Error(nVar + 1, ii, jj) = 0.
    END DO
END DO

IF (OutputFlag .EQ. 1) THEN
    CALL WriteDataToDisk_OCTAVE( nVar,(/nXs,nYs/),VarNameVisu,U_NVisu,tOut,iter,TRIM(FileNameOctU)//".dat",FVDim)
    !CALL WriteDataToDisk_OCTAVE( nVar,(/nXs,nYs/),VarNameVisu,Error,  tOut,iter,TRIM(FileNameOctE)//".dat",FVDim)
        
ELSEIF (OutputFlag .EQ. 2) THEN
    CALL WriteDataToDisk_TECPLOT(nVar,(/nXs,nYs/),MeshBary,VarNameVisu,U_NVisu,tOut,iter,TRIM(FileNameTecU)//".dat",FVDim)
    !CALL WriteDataToDisk_TECPLOT(nVar,(/nXs,nYs/),MeshBary,VarNameVisu,Error,  tOut,iter,TRIM(FileNameTecE)//".dat",FVDim)
    
ELSEIF (OutputFlag .EQ. 3) THEN
    CALL WriteDataToDisk_OCTAVE( nVar,(/nXs,nYs/),VarNameVisu,U_NVisu,tOut,iter,TRIM(FileNameOctU)//"_O.dat",FVDim)
    !CALL WriteDataToDisk_OCTAVE( nVar,(/nXs,nYs/),VarNameVisu,Error,  tOut,iter,TRIM(FileNameOctE)//"_O.dat",FVDim)
    CALL WriteDataToDisk_TECPLOT(nVar,(/nXs,nYs/),MeshBary,VarNameVisu,U_NVisu,tOut,iter,TRIM(FileNameTecU)//"_T.dat",FVDim)
    !CALL WriteDataToDisk_TECPLOT(nVar,(/nXs,nYs/),MeshBary,VarNameVisu,Error,  tOut,iter,TRIM(FileNameTecE)//"_T.dat",FVDim)    

ELSEIF (OutputFlag .EQ. 4) THEN
    CALL WriteDataToDisk_New(nVar,(/nXs,nYs/),MeshBary,U_NVisu,iter,TRIM(FileNameNewU)//".dat")
    !CALL WriteDataToDisk_New(nVar,(/nXs,nYs/),MeshBary,Error,  iter,TRIM(FileNameNewE)//".dat")    
ENDIF


!-------------------------------------------------------------------------------!
END SUBROUTINE WriteSolutionToDisk
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WriteMeshToDisk()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nXs, nYs, MeshBary, UNIT_FILE
USE MOD_FiniteVolume_vars,ONLY: Examplefolder
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
CHARACTER(LEN=255) :: FileName
!-------------------------------------------------------------------------------!
Filename = TRIM(Examplefolder//"/"//"MESH_CoordinateX.dat")
PRINT*, 'Writing mesh to disk: done'

OPEN(UNIT_FILE, FILE = Filename)
WRITE(UNIT_FILE,'(A)') "CoordinateX"

DO ii = 1,nXs
    WRITE(UNIT_FILE,'(SP,1(ES21.14E2))') MeshBary(1,ii,1)
END DO
CLOSE(UNIT_FILE)
!-------------------------------------------------------------------------------!
#ifdef TWODIM
Filename = TRIM(Examplefolder//"/"//"MESH_CoordinateY.dat")
WRITE(*,'(A10,1X,A)') " ", "Octave ASCII"//" => "//TRIM(FileName)
OPEN(UNIT_FILE, FILE = Filename)
WRITE(UNIT_FILE,'(A)') "CoordinateY"
DO jj = 1,nYs
    WRITE(UNIT_FILE,'(SP,1(ES21.14E2))') MeshBary(2,1,jj)
END DO
CLOSE(UNIT_FILE)
#endif
!-------------------------------------------------------------------------------!
END SUBROUTINE WriteMeshToDisk
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WriteDataToDisk_TECPLOT(nVar,nElems, Coords, VarNames,  &
                                   OutputData, tOut, iter, FileName, ZoneTitel)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nDims, UNIT_FILE
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN)          :: nVar
INTEGER,INTENT(IN)          :: nElems(2)
INTEGER,INTENT(IN)          :: iter
REAL,INTENT(IN)             :: tOut
REAL,INTENT(IN)             :: Coords(1:nDims,1:nElems(1),1:nElems(2))
REAL,INTENT(IN)             :: OutputData(1:nVar+1,1:nElems(1),1:nElems(2))
CHARACTER(LEN=*),INTENT(IN) :: VarNames(1:nVar+1)
CHARACTER(LEN=*),INTENT(IN) :: FileName
CHARACTER(LEN=*),INTENT(IN) :: ZoneTitel
!-------------------------------------------------------------------------------!
INTEGER                     :: ii, jj, iVar, Offset
CHARACTER(LEN=35)           :: VarString
CHARACTER(LEN=255)          :: FormatString
CHARACTER(LEN=512)          :: FormatTitle
CHARACTER(LEN=255)          :: temp1, temp2, temp3
REAL, DIMENSION(nVar+1)     :: tempOut
!-------------------------------------------------------------------------------!
PRINT*, achar(27)//'[92m'//REPEAT(' ',4)//'(Tecplot) '//TRIM(FileName)//achar(27)//'[0m'
FormatTitle = ""
#ifdef TWODIM
    Offset = 38
    FormatTitle(1:37) = 'VARIABLES="CoordinateX","CoordinateY"'
#else
    Offset = 24
    FormatTitle(1:23) = 'VARIABLES="CoordinateX"'
#endif
DO iVar = 1,nVar + 1
    !--------------------------------------------------------------------
    !Add the name of the variables to the FormatTitle string
    WRITE(VarString,'(A2,A,A1)') ',"', TRIM(VarNames(iVar)), '"'
    FormatTitle(Offset:Offset + LEN(TRIM(VarString))) = TRIM(VarString)
    !--------------------------------------------------------------------
    !Offset is to compute the length of "FormatTitle"
    Offset = Offset + LEN(TRIM(VarString))
    !--------------------------------------------------------------------
END DO
WRITE(temp1,'(F15.8)') tOut
WRITE(temp2,'(I8)')    nElems(1)
WRITE(temp3,'(I8)')    nElems(2)
!--------------------------------------------------------------------
OPEN(UNIT_FILE, FILE = TRIM(FileName), STATUS = "REPLACE")
WRITE(UNIT_FILE,'(A)') FormatTitle(1:Offset-1)
WRITE(UNIT_FILE,'(A,A,A,A,A,A,A,A)')                                           &
      'ZONE T="', TRIM(ZoneTitel), '",STRANDID=1, SOLUTIONTIME=',              &
      TRIM(ADJUSTL(TRIM(temp1))), ', DATAPACKING=POINT, ZONETYPE=ORDERED, I=', &
      TRIM(ADJUSTL(TRIM(temp2))), ', J = ', TRIM(ADJUSTL(TRIM(temp3)))
WRITE(FormatString,'(A4,I2,A15)') "(SP,", nDims + nVar + 1, "(ES21.14E2,2X))"
DO jj = 1, nElems(2)
    DO ii = 1, nElems(1)
        DO iVar = 1, nVar + 1
            IF (ABS(OutputData(iVar, ii, jj)) < 1.e-70) THEN
                tempOut(iVar) = 0.d0
            ELSE
                tempOut(iVar) = OutputData(iVar, ii, jj)
            END IF
        END DO
        WRITE (UNIT_FILE, FormatString) Coords(1:nDims, ii, jj), tempOut(1:nVar + 1)
        

    END DO
END DO
CLOSE(UNIT_FILE)
!-------------------------------------------------------------------------------!
END SUBROUTINE WriteDataToDisk_TECPLOT
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WriteDataToDisk_OCTAVE(nVar, nElems, VarNames,&
                                  OutputData, tOut, iter, FileName, ZoneTitel)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nDims, UNIT_FILE
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN)          :: nVar
INTEGER,INTENT(IN)          :: nElems(2)
REAL,INTENT(IN)             :: tOut
INTEGER,INTENT(IN)          :: iter
REAL,INTENT(IN)             :: OutputData(1:nVar+1,1:nElems(1),1:nElems(2))
CHARACTER(LEN=*),INTENT(IN) :: VarNames(1:nVar+1)
CHARACTER(LEN=*),INTENT(IN) :: FileName
CHARACTER(LEN=*),INTENT(IN) :: ZoneTitel
!-------------------------------------------------------------------------------!
INTEGER                     :: ii, jj, iVar
INTEGER                     :: Width, Pos
CHARACTER(LEN=35)           :: VarString
CHARACTER(LEN=255)          :: FormatString
CHARACTER(LEN=512)          :: FormatTitle
REAL, DIMENSION(nVar+1)     :: tempOut
!-------------------------------------------------------------------------------!


PRINT*, achar(27)//'[92m'//REPEAT(' ',4)//'(Octave)  '//TRIM(FileName)//achar(27)//'[0m'
Pos         = 1
Width       = 23
FormatTitle = ""
DO iVar = 1,nVar + 1
    WRITE(VarString,'(A)') TRIM(VarNames(iVar))
    FormatTitle(Pos:Pos + Width - 1) = ADJUSTL(TRIM(VarString))
    Pos = Pos + Width
END DO
OPEN(UNIT_FILE, FILE = TRIM(FileName), STATUS = "REPLACE")
WRITE(UNIT_FILE,'(A)') FormatTitle(1:Pos-1)
WRITE(FormatString,'(A4,I2,A15)') "(SP,", nVar + 1, "(ES21.14E2,2X))"
DO jj = 1, nElems(2)
    DO ii = 1, nElems(1)
        DO iVar = 1, nVar + 1
            IF (ABS(OutputData(iVar, ii, jj)) < 1.e-70) THEN
                tempOut(iVar) = 0.d0
            ELSE
                tempOut(iVar) = OutputData(iVar, ii, jj)
            END IF
        END DO
        WRITE (UNIT_FILE, FormatString) tempOut
    END DO
END DO
CLOSE(UNIT_FILE)
!-------------------------------------------------------------------------------!
END SUBROUTINE WriteDataToDisk_OCTAVE
!===============================================================================!
!
!
!
!
!===============================================================================!
SUBROUTINE WriteDataToDisk_NEW(nVar,nElems,Coords,OutputData,iter,FileName)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: nDims, UNIT_FILE
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! FORMAL ARGUMENTS
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN)          :: nVar
INTEGER,INTENT(IN)          :: nElems(2)
REAL,INTENT(IN)             :: Coords(1:nDims,1:nElems(1),1:nElems(2))
REAL,INTENT(IN)             :: OutputData(1:nVar+1,1:nElems(1),1:nElems(2))
INTEGER,INTENT(IN)          :: iter
CHARACTER(LEN=*),INTENT(IN) :: FileName
!-------------------------------------------------------------------------------!
INTEGER                     :: ii, jj, iVar
CHARACTER(LEN=255)          :: FormatString
REAL, DIMENSION(nVar+1)     :: tempOut
!-------------------------------------------------------------------------------!
ii = 0*iter

PRINT*, achar(27)//'[92m'//REPEAT(' ',4)//TRIM(FileName)//achar(27)//'[0m'
OPEN(UNIT_FILE, FILE = TRIM(FileName), STATUS = "REPLACE")
WRITE(FormatString,'(A4,I2,A15)') "(SP,", nDims + nVar + 1, "(ES21.14E2,2X))"
DO jj = 1, nElems(2)
    DO ii = 1, nElems(1)
        DO iVar = 1, nVar + 1
            IF (ABS(OutputData(iVar, ii, jj)) < 1.e-70) THEN
                tempOut(iVar) = 0.d0
            ELSE
                tempOut(iVar) = OutputData(iVar, ii, jj)
            END IF
        END DO
        WRITE (UNIT_FILE, FormatString) Coords(1:nDims, ii, jj), tempOut(1:nVar + 1)
    END DO
END DO
CLOSE(UNIT_FILE)
!-------------------------------------------------------------------------------!
END SUBROUTINE WriteDataToDisk_NEW
!===============================================================================!
! 
! 
! 
!===============================================================================!
SUBROUTINE ErrorCalculation
USE MOD_PhysicsFrame,ONLY: ExactFunction
USE MOD_FiniteVolume_vars
IMPLICIT NONE
!-------------------------------------------------------------------------------!
CHARACTER(LEN=255)                  :: FNAME
REAL, DIMENSION(nVar, nGPsX, nGPsY) :: Utemp
REAL, ALLOCATABLE                   :: Uexact(:,:,:), L1error(:)
INTEGER                             :: ii, jj, iGP, jGP, kk
!-------------------------------------------------------------------------------!
ALLOCATE( L1error(1:nVar) )
ALLOCATE( Uexact( 1:nVar, 1:nXs, 1:nYs))
!-------------------------------------------------------------------------------!

Uexact  = 0.
L1error = 0.
Utemp   = 0.

DO jj = 1, nYs
    DO ii = 1, nXs
        DO iGP = 1, nGPsX
            DO jGP = 1, nGPsY
                CALL ExactFunction(tEnd, MeshGP(:, ii, jj, iGP, jGP), Utemp(:, iGP, jGP))
                Uexact(:, ii, jj) = Uexact(:, ii, jj) + WeightsGP(iGP, jGP)*Utemp(:, iGP, jGP)
            END DO
        END DO
        L1error = L1error + ABS(U(:, ii, jj) - Uexact(:, ii, jj))
    END DO
END DO

DO kk = 1,nDims
   L1error = L1error * MESH_DX(kk)
END DO

FNAME = 'errorL1_tend.dat'
OPEN(666,FILE = Examplefolder//"/Error/"//TRIM(FNAME))
WRITE(666,*) L1error(1:nVar)
CLOSE(666)

DEALLOCATE( L1error )
DEALLOCATE( Uexact )

END SUBROUTINE ErrorCalculation
!===============================================================================!


END MODULE MOD_Output
!-------------------------------------------------------------------------------!
