!===============================================================================!
PROGRAM FiniteVolumeSWM
!-------------------------------------------------------------------------------!
USE MOD_Mesh,          ONLY: BuildMesh
USE MOD_Parameters,    ONLY: SimulationSelector, InitializeParameters
USE MOD_FiniteVolume,  ONLY: FillInitialConditions,  InitializeWBVariables
USE MOD_FiniteVolume,  ONLY: InitializeFiniteVolume, FinalizeFiniteVolume
USE MOD_Output,        ONLY: WriteMeshToDisk, WriteSolutionToDisk
USE MOD_Output,        ONLY: ErrorCalculation
!-------------------------------------------------------------------------------!
! For the time discretization:
USE MOD_TimeDiscretization,ONLY: TimeApproximation, TimeStep
USE MOD_FiniteVolume_vars, ONLY: dt, tEnd, dt_save, tGlobal
USE MOD_FiniteVolume_vars, ONLY: nOutputFiles, Examplefolder
USE MOD_FiniteVolume_vars, ONLY: nVar, nDims, nMoms, nGPsX, nXs, nYs, nGhostsX,    Identity, MESH_DX
USE MOD_FiniteVolume_vars, ONLY: OutputFlag, BathymetryFlag
USE MOD_FiniteVolume_vars, ONLY: SplittingFlag, ParameterNu, Hmin
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars, ONLY: VarTimeOpt, VarNConOpt, VarRecoOpt 
USE MOD_FiniteVolume_vars, ONLY: VarOutsOpt, VarInitOpt, VarBConOpt, VarModlOpt
USE MOD_FiniteVolume_vars, ONLY: VarBathOpt, VarFricOpt, auxBC, IndexPass, SSFriction

USE MOD_FiniteVolume_vars, ONLY: U,V
USE MOD_PhysicsFrame,      ONLY: PrimToCons, ConsToPrim
!-------------------------------------------------------------------------------!
USE MOD_Output,      ONLY: ReducedModelPostProcessing
USE NonlinearSolver, ONLY: NewtonRaphsonSolver
USE MOD_MomentModels

IMPLICIT NONE
!-------------------------------------------------------------------------------!
REAL    :: t, tsave, dt_min, inf_dt, sup_dt, start_time, end_time, elapsed_time,  Hminval
INTEGER :: iter, tot_iter, ii, jj, kk, ii0
REAL, ALLOCATABLE :: Utemp(:)
! testing
REAL, ALLOCATABLE :: Uaux(:,:,:), Vaux(:,:,:), Aplus(:,:), Amins(:,:)
REAL              :: dtdx, dxdt
!-------------------------------------------------------------------------------!
CALL InitializeParameters()
CALL SimulationSelector()
CALL InitializeFiniteVolume()

CALL BuildMesh()

iter = 0; t = 0.0;
dt_save = tEnd/REAL(nOutputFiles)
tsave   = t + dt_save
IF (t .EQ. tEnd) THEN
    RETURN
END IF
#ifdef WELLBALANCED
    CALL InitializeWBVariables()
#endif
CALL FillInitialConditions()

!-------------------------------------------------------------------------------!
IF (OutputFlag .EQ. 1 .OR. OutputFlag .EQ. 3) THEN
    CALL WriteMeshToDisk()
END IF

PRINT*, "Writing data to disk:" 
CALL WriteSolutionToDisk(t,iter)


OPEN(1, FILE = TRIM(Examplefolder)//"/timepoints.dat", STATUS = "REPLACE")

inf_dt   = 1.0E+10
sup_dt   = 0.0
tot_iter = 0
CALL CPU_TIME(start_time)

Hminval = 1.0E+10

ALLOCATE(Utemp(1:nVar))
Utemp = 0.0

ALLOCATE(Uaux(1:nVar,-1:nXs+1,1:nYs), Vaux(1:nVar,-1:nXs+1,1:nYs), Aplus(1:nVar,1:nVar), Amins(1:nVar,1:nVar))
Uaux  = 0.0
Vaux  = 0.0
Aplus = 0.0
Amins = 0.0

DO WHILE (t .LT. tEnd - 1.0E-10)
    
    tot_iter = tot_iter + 1
    CALL TimeStep(dt_min)
    !dt_min = 0.001 !! To fix delta t.
        
    dt     = MIN( MIN(dt_min,  tsave - t), MIN(dt_min,  tEnd - t) )
    inf_dt = MIN(inf_dt,dt)
    sup_dt = MAX(sup_dt,dt)

    CALL TimeApproximation(t)    

    IF ((SplittingFlag .EQ. 1)) THEN!-----------------------------------!  
        DO jj = 1, nYs
            DO ii = 1, nXs
                IndexPass = (/ii, jj/)
                IF (U(1,ii,jj) .LT. Hmin) THEN                
                    U(2:nVar,ii,jj) = 0.0   
                    V(2:nVar,ii,jj) = 0.0             
                ELSE 
                    ! Solve: (h^{n+1}I - dt*matrixS)*Prim^{n+1} = Cons^n
                    CALL NewtonRaphsonSolver(U(:,ii,jj), Utemp)
                    U(:,ii,jj) = Utemp
                    CALL ConsToPrim(U(:,ii,jj), V(:,ii,jj)) 
                    IF (U(1,ii,jj) .LT. Hmin) THEN                
                        U(2:nVar,ii,jj) = 0.0   
                        !V(2:nVar,ii,jj) = 0.0   
                    END IF
                END IF
            END DO
        END DO!----------------------------------------------------------!        
    ELSE IF (SplittingFlag .EQ. 0) THEN    
        DO jj = 1, nYs
            DO ii = 1, nXs
                IndexPass = (/ii, jj/)                
                IF (U(1,ii,jj) .LT. Hmin) THEN                
                    U(2:nVar,ii,jj) = 0.0   
                    V(2:nVar,ii,jj) = 0.0             
                ELSE                 
                    CALL SourceTerm(V(:,ii,jj), SSFriction(:,ii,jj))
                    U(:,ii,jj) = U(:,ii,jj) + dt*SSFriction(:,ii,jj)  
                    CALL ConsToPrim(U(:,ii,jj), V(:,ii,jj))                    
                    IF (U(1,ii,jj) .LT. Hmin) THEN                
                        U(2:nVar,ii,jj) = 0.0   
                        !V(2:nVar,ii,jj) = 0.0   
                    END IF
                END IF
            END DO
        END DO
    END IF!----------------------------------------------------------!
    
    
    t = t + dt
    tGlobal = t    
    
    ! To write the solution at the time tsave        
    IF (ABS(tsave - t) .LT. 1.0E-10) THEN
        iter = iter + 1             
        CALL WriteSolutionToDisk(t,iter)
        
        tsave = tsave + dt_save
        WRITE(1,*) t
        
        IF (tsave .GT. tEnd) THEN
            tsave = tEnd
        END IF
        
        IF (VarModlOpt .EQ. "ReducedModelI") THEN    
            CALL ReducedModelPostProcessing(iter)
        END IF
    
    END IF
    
END DO
CALL CPU_TIME(end_time)
elapsed_time = end_time - start_time
CLOSE(1)
!-------------------------------------------------------------------------------!
!CALL ErrorCalculation()  ! L1 error calculation at the end of the iteration.
CALL FinalizeFiniteVolume()

DEALLOCATE(Utemp)
DEALLOCATE(Uaux, Vaux, Aplus, Amins)

! Formal end of the numerical simulation -- The rest is the related information  
!-------------------------------------------------------------------------------!

IF (BathymetryFlag .GT. 0  ) THEN; VarBathOpt = "Yes"; ELSE; VarBathOpt = "No";  END IF
IF (ParameterNu    .EQ. 0.0) THEN; VarFricOpt = "No";  ELSE; VarFricOpt = "Yes"; ENDIF

!-------------------------------------------------------------------------------!
! FILE WITH INFORMATION:
OPEN(33, FILE = TRIM(Examplefolder)//"/information.txt", STATUS = "REPLACE")
WRITE(33,'(A28)')       "=========================== "
WRITE(33,'(A28)')       "INFORMATION FILE:         | " 
WRITE(33,'(A)')         "==========================|------------------------------------------------ "
WRITE(33,'(A28,I0,A)')  "Spatial dimension:          ", nDims, "D"
WRITE(33,'(A28,A)')     "Moment model:               ", TRIM(VarModlOpt)
WRITE(33,'(A28,I0)')    "Number of moments:          ", nMoms
WRITE(33,'(A28,A)')     "Friction (yes or no):       ", TRIM(VarFricOpt)
WRITE(33,'(A)')         "==========================|------------------------------------------------ "
WRITE(33,'(A28,A)')     "Time approximation:         ", TRIM(VarTimeOpt)
WRITE(33,'(A28,A)')     "Non-conservative approx.:   ", TRIM(VarNConOpt)
WRITE(33,'(A28,A)')     "-- Path used -----------:   ", "Linear"    
WRITE(33,'(A28,A)')     "Flux approx. (if any):      ", "Not implemented for now"
WRITE(33,'(A28,A)')     "Polynomial reconstruction:  ", TRIM(VarRecoOpt)
WRITE(33,'(A28,I0)')    "Number of Gauss Quad. P.:   ", nGPsX
WRITE(33,'(A28,I0)')    "Number of Ghost cells:      ", nGhostsX
WRITE(33,'(A28,A)')     "Bathymetry (yes or no):     ", TRIM(VarBathOpt)
WRITE(33,'(A28,A)')     "Initial Condition:          ", TRIM(VarInitOpt)

IF (nDims .EQ. 1) THEN
    WRITE(auxBC, "(A,A,A,A)") "left: ", TRIM(VarBConOpt(1)), &
                           ", right: ", TRIM(VarBConOpt(2))
ELSE
    WRITE(auxBC, "(A,A,A,A,A,A,A,A)") "left: ", TRIM(VarBConOpt(1)), &
                                   ", right: ", TRIM(VarBConOpt(2)), &
                                  ", bottom: ", TRIM(VarBConOpt(3)), &
                                     ", top: ", TRIM(VarBConOpt(4))                                       
END IF

WRITE(33,'(A28,A)')       "Boundary conditions:        ", TRIM(auxBC)
WRITE(33,'(A28,A)')       "Output type:                ", TRIM(VarOutsOpt)
WRITE(33,'(A)')           "==========================|------------------------------------------------ "
WRITE(33,'(A28,F10.8)')   "Minimum dt (s):             ", inf_dt
WRITE(33,'(A28,F10.8)')   "Maximum dt (s):             ", sup_dt
WRITE(33,'(A28,F10.8)')   "Simulation time (s):        ", elapsed_time
WRITE(33,'(A28,I0)')      "Total iterations:           ", tot_iter
WRITE(33,'(A28,I0)')      "Snapshots saved:            ", iter
WRITE(33,'(A28,I0)')      "Cells X (Cartesian):        ", nXs
WRITE(33,'(A28,I0)')      "Cells Y (Cartesian):        ", nYs
WRITE(33,'(A)')           "==========================|------------------------------------------------ "
WRITE(33,*) "------------------------------------------------------------------------------------------------------------------"
WRITE(33,*) '* Output dat files:                                                                                         '
WRITE(33,*) '    > The names starting with "soln" correspond to the saved solution;                                  '  
WRITE(33,*) '      it follows an integer, which is a counter of the saved snapshots (not the time iteration);            '
WRITE(33,*) '      the spatial dimension is placed at the end preceded by an underscore.\n                               '
SELECT CASE (OutputFlag)
CASE(1)
    WRITE(33,*) '    > The format of the solution files is made to be read by Octabe or matlab through the "load" function:'
    WRITE(33,*) '      - The first row contains format information'
    WRITE(33,'(A15,I0,A)') '       - First ', nVar,' columns: Each component of the solution,'
CASE(2)
    WRITE(33,*) '    > The format of the solution files, is made to be read by the software Tecplot:'
    WRITE(33,*) '      - The first and second row contain format information'
CASE(3)
    WRITE(33,*) '    > The format of both solution files are made to be read by software Tecplot and as a matrix, both follow the format:'
    WRITE(33,*) '      - Octave files "_O":  The first row contains format information'
    WRITE(33,'(A15,I0,A)') '       - First ', nVar,' columns: Each component of the solution,'
    WRITE(33,*)            '      - Last column:  Bathymetry or Topograph                                                       '
    WRITE(33,*)            '      ---------------------'
    WRITE(33,*) '      - Tecplot files "_T": The first and second row contain format information'
CASE(4)
    WRITE(33,*) '    > The format of the solution files is made to be read as a matrix'
END SELECT
!----------------------------------------------------------------------------------------------------------------
IF ((nDims .EQ. 1 ).AND. (OutputFlag .GT. 1)) THEN   
    WRITE(33,*)            '      - First column: x-coordinate'  
    WRITE(33,'(A19,I0,A)') '       - Second to ',nVar+1,' columns: Each component of the solution'
ELSE IF ((nDims .EQ. 2 ).AND. (OutputFlag .GT. 1))  THEN
        WRITE(33,*)            '      - First two columns: x-coordinate and y-coordinate'  
        WRITE(33,'(A18,I0,A)') '       - Third to ',nVar+2,' columns: Each component of the solution'    
END IF
!----------------------------------------------------------------------------------------------------------------
WRITE(33,*) '      - Last column:  Bathymetry or Topograph                                                       '
WRITE(33,*) '      - Each row: Evaluation of the variables at each node of the mesh.                             '
SELECT CASE (OutputFlag)
CASE(1,3)
    WRITE(33,*) '    \n> Files corresponding to the mesh coordinates are generated for to run the octave files.      '
END SELECT
WRITE(33,*) "------------------------------------------------------------------------------------------------------------------"
WRITE(33,*) '* The file "timepoints.dat" contains the time points in which the solutions files have been saved,'
WRITE(33,*) '  the order of these times is in agreement with the number of the saved file.                     '
WRITE(33,*) "------------------------------------------------------------------------------------------------------------------"
WRITE(33,*) '* The produced error files saved in the Error folder contain the corresponding error of the saved solution'
WRITE(33,*) '  with respect to a given exact solution.'
WRITE(33,*) "------------------------------------------------------------------------------------------------------------------"
        
CLOSE(33)
PRINT*, REPEAT('-',50)
PRINT *, achar(27)//'[1;42m'//' Finish. Find all details in information.dat.     '//achar(27)//'[0m'
PRINT*, REPEAT('-',50)
!-------------------------------------------------------------------------------!
END PROGRAM FiniteVolumeSWM
!===============================================================================!


