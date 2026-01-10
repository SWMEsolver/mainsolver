!===============================================================================!
MODULE MOD_Mesh
!===============================================================================!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE BuildMesh
    MODULE PROCEDURE BuildMesh
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: BuildMesh
!-------------------------------------------------------------------------------!
CONTAINS
!===============================================================================!


!===============================================================================!
SUBROUTINE BuildMesh()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY: MESH_SX, MESH_X0, MESH_X1, MESH_DX
USE MOD_FiniteVolume_vars,ONLY: nDims, nGPsX, nGPsY, nXs, nYs
USE MOD_FiniteVolume_vars,ONLY: NormVectX, NormVectY, TangVectX, TangVectY
USE MOD_FiniteVolume_vars,ONLY: MeshNodes, MeshBary, MeshGP
USE MOD_FiniteVolume_vars,ONLY: WeightsGP, WeightsGPBnd    
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
INTEGER :: ii, jj, iGP, jGP
REAL, DIMENSION(nGPsX) :: quadWeights1D, quadNodes1D 
REAL :: auxnxy(1:2)
REAL :: auxidx(1:2)
REAL :: auxGPvec(1:2)

auxnxy = (/REAL(nXs), REAL(nYs)/)
MeshNodes = 0.0
MeshBary  = 0.0
Mesh_SX   = ABS(Mesh_X1 - Mesh_X0)
DO ii = 1,nDims
    Mesh_DX(ii) = ABS(Mesh_SX(ii))/auxnxy(ii)
END DO

DO jj = 0, nYs
    DO ii = 0, nXs
        auxidx = (/REAL(ii), REAL(jj)/)
        MeshNodes(:, ii, jj) = Mesh_X0(1:nDims) + auxidx(1:nDims)*Mesh_DX
    END DO
END DO
DO jj = 1, nYs
    DO ii = 1, nXs
        MeshBary(:, ii, jj) = MeshNodes(:, ii - 1, jj - 1) + 0.5*Mesh_DX
    END DO
END DO

!---------------------------!
! Only used in the 2D case
DO iGP = 1, nGPsY
    DO jj = 1, nYs
        DO ii = 0, nXs
            NormVectX(:, iGP, ii, jj) = (/1.0, 0.0/)
            TangVectX(:, iGP, ii, jj) = (/0.0, 1.0/)
        END DO
    END DO
    DO jj = 0, nYs
        DO ii = 1, nXs
            NormVectY(:, iGP, ii, jj) = (/0.0, 1.0/)
            TangVectY(:, iGP, ii, jj) = (/-1.0, 0.0/)
        END DO
    END DO
END DO

!------------------------------!
!   MeshGP Quadrature          !
!------------------------------!
! These are classical 1D-Gauss quadrature weights and points.
! See https://engcourses-uofa.ca/books/numericalanalysis/numerical-integration/gauss-quadrature/
SELECT CASE (nGPsX)
CASE (1)
    quadWeights1D(1) = 1.0
    quadNodes1D(1)   = 0.0
CASE (2)
    quadWeights1D = (/0.5, 0.5/)
    quadNodes1D   = (/-1./(2.*sqrt(3.)), 1./(2.*sqrt(3.))/)
CASE (3)
    quadWeights1D = (/5./18., 4./9., 5./18./)
    quadNodes1D   = (/-0.5*sqrt(3./5.), 0., 0.5*sqrt(3./5.)/)
CASE (4)
    quadWeights1D = (/(18. - SQRT(30.))/72., (18. + SQRT(30.))/72., (18. + SQRT(30.))/72., (18. - SQRT(30.))/72./)
    quadNodes1D   = (/-0.5*SQRT(3./7. + 2./7.*SQRT(6./5.)), -0.5*SQRT(3./7. - 2./7.*SQRT(6./5.)), &
                       0.5*SQRT(3./7. - 2./7.*SQRT(6./5.)),  0.5*SQRT(3./7. + 2./7.*SQRT(6./5.))/)
CASE DEFAULT
    PRINT *, "Quadrature not implemented"
    STOP
END SELECT

DO iGP = 1, nGPsX
    WeightsGPBnd(iGP) = quadWeights1D(iGP)
    DO jGP = 1, nGPsY
        IF (nDims .EQ. 2) THEN
            WeightsGP(iGP, jGP) = quadWeights1D(iGP)*quadWeights1D(jGP)
        ELSE
            WeightsGP(iGP, jGP) = quadWeights1D(iGP) ! Only on the first argument, iGP
        END IF        
        DO jj = 1, nYs
            DO ii = 1, nXs            
                auxGPvec = (/ quadNodes1D(iGP) , quadNodes1D(jGP) /)
                MeshGP(1:nDims, ii, jj, iGP, jGP) = MeshBary(1:nDims, ii, jj) + auxGPvec(1:nDims)*Mesh_DX(1:nDims)                
            END DO
        END DO
    END DO
END DO

END SUBROUTINE BuildMesh

!===============================================================================!
END MODULE MOD_Mesh
!-------------------------------------------------------------------------------!
