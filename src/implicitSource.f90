!===============================================================================!
MODULE MOD_ImplicitSource
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ImplicitSource2
  MODULE PROCEDURE ImplicitSource2
END INTERFACE

!-------------------------------------------------------------------------------!
PUBLIC :: ImplicitSource2
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
SUBROUTINE ImplicitSource2(Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume_vars,ONLY:  Hmin

USE MOD_FiniteVolume_vars,ONLY: ParameterNu
USE MOD_FiniteVolume_vars,ONLY: ParameterLambda
USE MOD_FiniteVolume_vars,ONLY: dt

!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL :: CTE1, CTE2

REAL,INTENT(INOUT) :: Cons(1:4)

REAL             :: ht
REAL             :: InvA(1:4, 1:4)
REAL             :: res(1:4)
INTEGER          :: ii,jj,kk
!REAL             :: detC
integer, parameter :: n = 4
!integer :: info, lda, ldvl, ldvr, lwork
!double precision :: wr(n), wi(n)
!double precision, allocatable :: work(:)
!double precision :: vl(n,n), vr(n,n)  ! Left and right eigenvectors

!-------------------------------------------------------------------------------!
CTE1 = ParameterNu
CTE2 = ParameterLambda
ht   = Cons(1) 
ii = 0
jj = 0
kk = 0
InvA = 0.0
res = 0.0


InvA(1,1) = 1.

InvA(2,2) = (CTE2*ht**5. + 240.*CTE1**2.*dt**2.*ht*(3.*CTE2 + ht) + 8.*CTE1*dt*ht**3.*(9.*CTE2 + ht))/&
     -  (720.*CTE1**3.*dt**3. + CTE2*ht**5. + 9.*CTE1*dt*ht**3.*(8.*CTE2 + ht) + 24.*CTE1**2.*dt**2.*ht*(30.*CTE2 + 13.*ht))


InvA(2,3) = -((CTE1*dt*ht**2.*(60.*CTE1*dt + ht**2.))/(720.*CTE1**3.*dt**3. + CTE2*ht**5. + 9.*CTE1*dt*ht**3.*(8.*CTE2 + ht) + 24.*CTE1**2.*dt**2.*ht*(30.*CTE2 + 13.*ht)))
InvA(2,4) = -((CTE1*dt*ht**2.*(12.*CTE1*dt + ht**2.))/(720.*CTE1**3.*dt**3. + CTE2*ht**5. + 9.*CTE1*dt*ht**3.*(8.*CTE2 + ht) + 24.*CTE1**2.*dt**2.*ht*(30.*CTE2 + 13.*ht)))


InvA(3,2) = (-3.*CTE1*dt*ht**2.*(60.*CTE1*dt + ht**2.))/(720.*CTE1**3.*dt**3. + CTE2*ht**5. + 9.*CTE1*dt*ht**3.*(8.*CTE2 + ht) + 24.*CTE1**2.*dt**2.*ht*(30.*CTE2 + 13.*ht))
InvA(3,3) =    (ht**2.*(60.*CTE1**2.*dt**2. + CTE2*ht**3. + 6.*CTE1*dt*ht*(10.*CTE2 + ht)))/&
     -  (720.*CTE1**3.*dt**3. + CTE2*ht**5. + 9.*CTE1*dt*ht**3.*(8.*CTE2 + ht) + 24.*CTE1**2.*dt**2.*ht*(30.*CTE2 + 13.*ht))
InvA(3,4) = (-3.*CTE1*dt*ht**4.)/(720.*CTE1**3.*dt**3. + CTE2*ht**5. + 9.*CTE1*dt*ht**3.*(8.*CTE2 + ht) + 24.*CTE1**2.*dt**2.*ht*(30.*CTE2 + 13.*ht))

InvA(4,2) = (-5.*CTE1*dt*ht**2.*(12.*CTE1*dt + ht**2.))/(720.*CTE1**3.*dt**3. + CTE2*ht**5. + 9.*CTE1*dt*ht**3.*(8.*CTE2 + ht) + 24.*CTE1**2.*dt**2.*ht*(30.*CTE2 + 13.*ht))
InvA(4,3) = (-5.*CTE1*dt*ht**4.)/(720.*CTE1**3.*dt**3. + CTE2*ht**5. + 9.*CTE1*dt*ht**3.*(8.*CTE2 + ht) + 24.*CTE1**2.*dt**2.*ht*(30.*CTE2 + 13.*ht))
InvA(4,4) = (ht**2.*(12.*CTE1**2.*dt**2. + CTE2*ht**3. + 4.*CTE1*dt*ht*(3.*CTE2 + ht)))/&
     -  (720.*CTE1**3.*dt**3. + CTE2*ht**5. + 9.*CTE1*dt*ht**3.*(8.*CTE2 + ht) + 24.*CTE1**2.*dt**2.*ht*(30.*CTE2 + 13.*ht))

!lda = n
!ldvl = n
!ldvr = n
!lwork = 4*n
!allocate(work(lwork))
! Compute eigenvalues only (WR, WI). VL and VR not referenced.
!call dgeev('N', 'N', n, invA, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

!IF (tot_iter .GT. 400) THEN
!if (info /= 0) then
!       PRINT *, "Error: DGEEV failed with info = ", info
!    else
!        PRINT *, "Eigenvalues:"
!        do ii = 1, n
!            if (wi(ii) == 0.0d0) then
!                PRINT*, "λ", ii, " = ", wr(ii)
!            else
!                PRINT*, "λ", ii, " = ", wr(ii), " + ", wi(ii), "ii"
!            end if
!        end do
!    end if
!END IF

!PRINT*, Cons
DO ii = 1, 4
    res(ii) = 0.0
    DO jj = 1, 4
      res(ii) = res(ii) + InvA(ii, jj) * Cons(jj)
  END DO
END DO

!PRINT*, "Result" 
!PRINT*, res(2)

DO kk = 1, 4
    Cons(kk)=res(kk)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE ImplicitSource2
!===============================================================================!
!
!
END MODULE MOD_ImplicitSource
!===============================================================================!
