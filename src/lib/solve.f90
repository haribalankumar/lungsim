
module solve
!*Brief Description:* This module contains solvers required for lung problems
!
!*LICENSE:*
!
!
!
!*Full Description:*
!
!Solvers included in this module are:
! BICGSTAB == The STABalized BIConjugate Gradient method
! GMRES == Generalised Minimal RESidual method
  use math_utilities
  use arrays,only: dp
  implicit none
  private
  public BICGSTAB_LinSolv
  public pmgmres_ilu_cr
  public solve_linear_system
contains
!
!#######################################################################
!
!*BICGSTAB_LinSolv:*
!   This function uses a preconditioned biconjugate gradient method to solve a
!   linear system Ax=b. Currently this solver is hard-coded to use the diaganol
!   entries of A as a preconditioner (the Jacobi preconditioner). However, this
!   could easily be added as an input variable. Sparsity is defined by compressed
!   row storage.
! Created by ARC: 15-12-2015

!   Input:
!   MatrixSize -> the size of the matrix.
!   NonZeros -> The number of non-zero entries.
!   SparseCol -> Define the column indices that have non-zero entries.
!   SparseRow ->Define the indices of SparseCol that represent each row.
!   SparseVal -> The values of non-zero entries.
!   RHS -> The vector on the RHS of the linear system

!   FLAGS - 0=Solution found to tolerance
!           1=No convergence max iterations reached
!           -1=breakdown rho=0
!           -2=breakdown omega=0

subroutine BICGSTAB_LinSolv(MatrixSize,NonZeros,RHS,Solution,SparseCol,SparseRow,SparseVal,TOLER,MaxIter)
    use arrays

    !Input/Output Variables
    integer, intent(in) :: MatrixSize,NonZeros
    integer, intent(in) :: SparseCol(NonZeros),SparseRow(MatrixSize+1),MaxIter
    double precision :: Solution(MatrixSize),SparseVal(NonZeros),RHS(MatrixSize),TOLER

    !Local Variables
    integer :: i,k,j,kstart,kend,Iter,FLAG
    double precision :: Ax(MatrixSize),bnrm2,error_r,PreCon(MatrixSize),PreConInv(MatrixSize)
    double precision :: r(MatrixSize),rtilde(MatrixSize),alpha,beta,omega,rho,rho_1
    double precision :: p(MatrixSize),v(MatrixSize),phat(MatrixSize),s(MatrixSize),snorm,resid,shat(MatrixSize),t(MatrixSize)
    real :: start,finish
    call CPU_TIME(start)
    print *, 'In Solver, Tolerance=', TOLER, 'Max iterations=', MaxIter

!!! PRECONDITIONER MATRIX - THE DIAGONALS OF THE MATRIX A P=diag(A), Jacobi preconditioner Inverse is deltaij/Aij
    DO i=1,MatrixSize
        kstart=SparseRow(i)
        kend=SparseRow(i+1)-1
        DO k=kstart,kend
            IF(SparseCol(k).eq.i)THEN
                PreCon(i)=SparseVal(k)
                PreConInv(i)=1/PreCon(i)
            ENDIF
        ENDDO
    ENDDO


    !bnrm2=norm(RHS)
    bnrm2=sqrt(DOT_PRODUCT(RHS,RHS))!Checked against matlab
    IF(bnrm2.eq.0.0_dp) bnrm2=1.0_dp


    error_r=0.0_dp
    !CALL ax_cr(MatrixSize,SparseRow,SparseCol,SparseVal,solution,Ax)
    DO i=1,MatrixSize !This loop depends on sparsity structure
        Ax(i)=0.0_dp !Initialise Ax for this entry
        kstart=SparseRow(i)
        kend=SparseRow(i+1)-1
        DO k=kstart,kend
            Ax(i)=Ax(i)+SparseVal(k)*solution(SparseCol(k)) !A*X
        ENDDO
    ENDDO
    r=RHS-Ax
    rtilde=r

    error_r=sqrt(DOT_PRODUCT(r,r))/bnrm2 !Checked against matlab


    IF(error_r.lt.TOLER)THEN
        write(*,*) 'Initial guess is perfect, Exiting solver',TOLER
        FLAG=0
    RETURN
    ENDIF
    alpha=0.0_dp
    beta=0.0_dp
    omega=1.0_dp

!!START OF ITERATIVE LOOP
     iterative_loop: DO Iter=1,MaxIter
        !rho is dot product of r and t tilde
        rho=DOT_PRODUCT(r,rtilde) !checked against matlab
        IF(rho.EQ.0.0_dp) EXIT iterative_loop
        IF(ITER.GT.1.0_dp)THEN
            beta=(rho/rho_1)*(alpha/omega)
            p=r+beta*(p-omega*v)
            phat=PreConInv*p
            v=0.0_dp*v
        ELSE
            p=r
            phat=PreConInv*p
            v=0.0_dp*p
        ENDIF


        !v=ApHat
        DO i=1,MatrixSize !This loop depends on sparsity structure
            kstart=SparseRow(i)
            kend=SparseRow(i+1)-1
            DO k=kstart,kend
                v(i)=v(i)+SparseVal(k)*phat(SparseCol(k))
            ENDDO
        ENDDO
        !v is correct against matlab

        alpha=rho/DOT_PRODUCT(rtilde,v) !alpha is correct first iter
        s=r-alpha*v

        snorm=sqrt(DOT_PRODUCT(s,s))!snorm correct against matlab

        !!EARLY CONVERGENCE TEST
        IF(snorm.LE.TOLER)THEN
            solution=solution+alpha*phat
            resid=snorm/bnrm2
            EXIT iterative_loop
        ENDIF


        !!Stabiliser shat=precon-1*s
        shat=PreConInv*s
        t=0.0_dp*s
        DO i=1,MatrixSize
            kstart=SparseRow(i)
            kend=SparseRow(i+1)-1
            DO k=kstart,kend
                t(i)=t(i)+SparseVal(k)*shat(SparseCol(k)) !A*X !correct first iter
            ENDDO
        ENDDO
        omega=DOT_PRODUCT(t,s)/DOT_PRODUCT(t,t)!omega corect first iter


        solution=solution+alpha*phat+omega*shat
        r=s-omega*t
        error_r=SQRT(DOT_PRODUCT(r,r))/bnrm2

        IF(error_r.LE.TOLER) EXIT iterative_loop
        IF(omega.EQ.0.0_dp) EXIT iterative_loop
        rho_1=rho!

    ENDDO iterative_loop
    IF(error_r.LE.TOLER.OR.snorm.LE.TOLER)THEN
        IF(snorm.LE.TOLER)THEN
            error_r=SQRT(DOT_PRODUCT(s,s))/bnrm2
        ENDIF
!        FLAG=0
        print *, 'Converged with error=',error_r
    ELSEIF(omega.EQ.0)THEN
!        FLAG=-2
        print *, 'Not converged omega =0'
    ELSEIF(rho.EQ.0)THEN
!        FLAG=-1
        print *, 'Not converged rho =0'
    ELSE
!        FLAG=1
        print *, 'Not converged max no of iterations reached'
        print *, 'Error=',error_r
    ENDIF
    CALL CPU_TIME(finish)
    print '("Time= ",f12.9,"seconds")',finish-start



end subroutine BICGSTAB_LinSolv
!
!#############################################################################
!
!*GMRES solver:*
!!! The following is the code from mgmres.f90, downloaded from
!!! www.people.sc.fsu.edu/~jbukhardt/f_src/mgmres/mgmres.html on 13/01/2016

!!!  Licensing:
    !    This code is distributed under the GNU LGPL license.
!!!  Modified:
    !    17 July 2007
!!!  Author:
    !    Original C version by Lili Ju.
    !    FORTRAN90 version by John Burkardt.
!!!  Reference:
    !    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    !    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    !    Charles Romine, Henk van der Vorst,
    !    Templates for the Solution of Linear Systems:
    !    Building Blocks for Iterative Methods,
    !    SIAM, 1994.
    !    ISBN: 0898714710,
    !    LC: QA297.8.T45.
    !
    !    Tim Kelley,
    !    Iterative Methods for Linear and Nonlinear Equations,
    !    SIAM, 2004,
    !    ISBN: 0898713528,
    !    LC: QA297.8.K45.
    !
    !    Yousef Saad,
    !    Iterative Methods for Sparse Linear Systems,
    !    Second Edition,
    !    SIAM, 2003,
    !    ISBN: 0898715342,
    !    LC: QA188.S17.

!!!  Parameters:
    !    Input, integer ( kind = 4 ) N, the order of the system.
    !
    !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
    !
    !    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
    !    indices of the matrix values.  The row vector has been compressed.
    !
    !    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
    !
    !    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A'.
    !
    !    Output, real ( kind = 8 ) W(N), the value of A'*X.

!!! Information:
    !    The Sparse Compressed Row storage format is used.
    !
    !    The matrix A is assumed to be sparse.  To save on storage, only
    !    the nonzero entries of A are stored.  The vector JA stores the
    !    column index of the nonzero value.  The nonzero values are sorted
    !    by row, and the compressed row vector IA then has the property that
    !    the entries in A and JA that correspond to row I occur in indices
    !    IA[I] through IA[I+1]-1.


  !*****************************************************************************




  !*****************************************************************************
  subroutine pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, &
       tol_abs, tol_rel,FLAG )
!!! PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
    !    This routine uses the incomplete LU decomposition for the
    !    preconditioning.  This preconditioner requires that the sparse
    !    matrix data structure supplies a storage position for each diagonal
    !    element of the matrix A, and that each diagonal element of the
    !    matrix A is not zero.
    !    Input/output, real ( kind = 8 ) X(N); on input, an approximation to
    !    the solution.  On output, an improved approximation.
    !
    !    Input, real ( kind = 8 ) RHS(N), the right hand side of the linear system.
    !
    !    Input, integer ( kind = 4 ) ITR_MAX, the maximum number of (outer)
    !    iterations to take.
    !
    !    Input, integer ( kind = 4 ) MR, the maximum number of (inner) iterations
    !    to take.  MR must be less than N.
    !
    !    Input, real ( kind = 8 ) TOL_ABS, an absolute tolerance applied to the
    !    current residual.
    !
    !    Input, real ( kind = 8 ) TOL_REL, a relative tolerance comparing the
    !    current residual to the initial residual.
    !
    implicit none

    integer ( kind = 4 ) mr
    integer ( kind = 4 ) n
    integer ( kind = 4 ) nz_num
    integer ( kind = 4 ) ia(*) !ia(n+1)
    integer ( kind = 4 ) ja(*) !ja(nz_num)
    integer ( kind = 4 ) itr_max

    real ( kind = 8 ) a(*) !a(nz_num)
    real ( kind = 8 ) x(*) !x(n)
    real ( kind = 8 ) rhs(*) !rhs(n)
    real ( kind = 8 ) tol_abs
    real ( kind = 8 ) tol_rel
    integer (kind = 4) FLAG


    real ( kind = 8 ) av
    real ( kind = 8 ) c(mr+1)
    real ( kind = 8 ), parameter :: delta = 1.0D-03
    real ( kind = 8 ) g(mr+1)
    real ( kind = 8 ) h(mr+1,mr)
    real ( kind = 8 ) htmp
    integer ( kind = 4 ) i
    integer ( kind = 4 ) itr
    integer ( kind = 4 ) itr_used
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) k_copy

    real ( kind = 8 ) l(ia(n+1)+1)
    real ( kind = 8 ) mu
    real ( kind = 8 ) r(n)
    real ( kind = 8 ) rho
    real ( kind = 8 ) rho_tol
    real ( kind = 8 ) s(mr+1)
    integer ( kind = 4 ) ua(n)
    real ( kind = 8 ) v(n,mr+1);
    logical, parameter :: verbose = .false.
    real ( kind = 8 ) y(mr+1)


    itr_used = 0
    FLAG=0 !not converged 

    call rearrange_cr ( n, ia, ja, a )

    call diagonal_pointer_cr ( n, ia, ja, ua )

    call ilu_cr ( n, nz_num, ia, ja, a, ua, l )

    if ( verbose ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'PMGMRES_ILU_CR'
       write ( *, '(a,i4)' ) '  Number of unknowns = ', n
    end if

    do itr = 1, itr_max

       call ax_cr ( n, ia, ja, a, x, r )

       r(1:n) = rhs(1:n) - r(1:n)

!!!! note that uses r,r in the call, but r,z in subroutine ??
       call lus_cr ( n, ia, ja, l, ua, r, r )

       rho = sqrt ( dot_product ( r, r ) )

       if ( verbose ) then
          write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
       end if

       if ( itr == 1 ) then
          rho_tol = rho * tol_rel
       end if

       v(1:n,1) = r(1:n) / rho

       g(1) = rho
       g(2:mr+1) = 0.0D+00

       h(1:mr+1,1:mr) = 0.0D+00

       do k = 1, mr

          k_copy = k

          call ax_cr ( n, ia, ja, a, v(1:n,k), v(1:n,k+1) )

          call lus_cr ( n, ia, ja, l, ua, v(1:n,k+1), v(1:n,k+1) )

          av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

          do j = 1, k
             h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
             v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
          end do

          h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

          if ( ( av + delta * h(k+1,k)) == av ) then
             do j = 1, k
                htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
                h(j,k) = h(j,k) + htmp
                v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
             end do
             h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
          end if

          if ( h(k+1,k) /= 0.0D+00 ) then
             v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
          end if

          if ( 1 < k ) then
             y(1:k+1) = h(1:k+1,k)
             do j = 1, k - 1
                call mult_givens ( c(j), s(j), j, y )
             end do
             h(1:k+1,k) = y(1:k+1)
          end if

          mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

          c(k) = h(k,k) / mu
          s(k) = -h(k+1,k) / mu
          h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
          h(k+1,k) = 0.0D+00
          call mult_givens ( c(k), s(k), k, g )

          rho = abs ( g(k+1) )

          itr_used = itr_used + 1

          if ( verbose ) then
             write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
          end if

          if ( rho <= rho_tol .and. rho <= tol_abs ) then
             FLAG=1
             exit
          elseif (isnan(rho))then
             FLAG=2
             exit
          end if

       end do

       k = k_copy - 1

       y(k+1) = g(k+1) / h(k+1,k+1)

       do i = k, 1, -1
          y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
       end do

       do i = 1, n
          x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
       end do

       if ( rho <= rho_tol .and. rho <= tol_abs ) then
          FLAG=1
          exit
       elseif (isnan(rho))then
          FLAG=2
          exit
       end if

    end do
    


    if ( verbose ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
       write ( *, '(a,i6)' ) '  Iterations = ', itr_used
       write ( *, '(a,g14.6)' ) '  Final residual = ', rho
    end if

    return
  end subroutine pmgmres_ilu_cr

  !*****************************************************************************


  subroutine solve_linear_system(LDA,M,N,num_depvar,A,B,X,FIRST_A)

!!! dummy arguments
    integer :: LDA,M,N,num_depvar
    double precision :: A(:),B(:),X(:)

    logical :: first_a
!!! local variables
    integer :: SOLVER !,ISOLV1(200000)
    integer,allocatable :: isolv1(:)
    double precision  :: SOLV_ANORM
    logical :: FACTOR
    
    allocate(isolv1(num_depvar*num_depvar))

    FACTOR=.TRUE.
    SOLVER = 1 ! for LAPACK
    !       Dense solvers: use LAPACK
    if(FACTOR) then
       CALL DIR_FACTOR(A,LDA,M,N,ISOLV1,SOLV_ANORM)
       FIRST_A=.FALSE.
    endif
    !         Solve the problem
    CALL DIR_SOLVE(A,LDA,LDA,M,N,X,B,ISOLV1)

    deallocate(isolv1)
    
  end subroutine solve_linear_system

  
!!! ##########################################################################  

subroutine DIR_FACTOR(A,LDA,M,N,IWORK_A,ANORM)
    
    !!!    DIR_FACTOR factorises a system that is to be used in
    !!!    one of the dense direct solvers
    !     Parameter List
    integer :: LDA,M,N,IWORK_A(*)
    double precision :: A(*),ANORM

    !     Local Constants
    integer NORMTYPE,SPARSE_A
    parameter(NORMTYPE=2,SPARSE_A=0)
    !     Local variables
    integer INFO,ISR_A(1)
    character(len=6) :: NAME
    
    
    INFO=0
    NAME=' '
    
    !     Get the norm of A
    CALL CALC_MAT_NORM(ANORM,A,LDA,N,SPARSE_A,ISR_A,NORMTYPE)
    
    !     Factorise the system:
    
    !     LU Decomposition
    CALL DGETRF(M,N,A,LDA,IWORK_A,INFO)
    NAME='DGETRF'
    
  end subroutine DIR_FACTOR
  
!!! ##########################################################################   

subroutine DIR_SOLVE(A,LDA,MMAX,M,N,X,B,IWORK_A)
    
    !!!    DIR_SOLV is a wrapper around the direct solvers, 
    !!!    allowing a more intelligable naming of arrays etc.
    !     Parameter List
    integer IWORK_A(*),LDA,MMAX,M,N
    double precision :: A(*),X(*),B(*)
    !     Local parameters
    integer INFO,NRHS,LDVT
    double precision :: SVDTOL
    character*(*) TRANS*(1),NAME*(6)
    !     Functions
    !double precision :: DNRM2
!joyce commented
    !EXTERNAL DCOPY,DNRM2
    
    
    TRANS='N'
    NRHS=1
    LDVT=N
    SVDTOL=1.0e-7_dp
    INFO=0
    NAME=' '
    
    !     LU Decomposition
    CALL DCOPY(M,B,1,X,1)
    CALL DGETRS(TRANS,N,NRHS,A,LDA,IWORK_A,X,MMAX,INFO)
    NAME='DGETRS'
    
    !     Catch errors from LAPACK routines
    if(INFO.NE.0) then
       WRITE(*,'('' >>INFO='',I6,'' in  '',A,'''')') INFO,NAME
    endif
    
  end subroutine DIR_SOLVE
  

!!! ##########################################################################   
     
  
  subroutine CALC_MAT_NORM(ANORM,A,LDA,N,SPARSE_A,ISR_A,NORMTYPE)
    
!!!    Calculates the infinity, first or second norms of a matrix,
!!!    ||A||_\infty, ||A||_1 or ||A||_2. Used in determining the
!!!    "scale" of a matrix, when considering convergance with 
!!!    iterative solver.

!!! dummy arguments
    integer LDA,N,SPARSE_A,ISR_A(:),NORMTYPE
    real(dp) :: ANORM,A(*)
!!! local variables
    integer :: I,IJ
    real(dp) :: TMP,SUM
!!! external functions (solver library)
    !real(dp) :: DASUM,PDASUM ,DNRM2
    !EXTERNAL DNRM2,DASUM,PDASUM
    
    !     Infinity norm
    if(NORMTYPE.EQ.0) then     !  Non-sparse matricies
       if(SPARSE_A.EQ.0) then
          TMP=0.0_dp
          do I=1,N
             TMP=MAX(TMP,DASUM(N,A(I),LDA))
          enddo
          ANORM=TMP
       else if(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) then   !Sparse matrices: CMISS compressed row
          TMP=0.0_dp
          do I=1,N
             SUM=0.0_dp
             do IJ=ISR_A(I),ISR_A(I+1)-1
                SUM=SUM+ABS(A(IJ))
             enddo
             TMP=MAX(TMP,SUM)
          enddo
          ANORM=TMP
          
       else
          !          ERROR='>>Sparsity type not implemented'
          !          GOTO 9999
       endif
       
       !     Second norm
    else if(NORMTYPE.EQ.2) then
       
       if(SPARSE_A.EQ.0) then  !       Non-sparse matrices
          SUM=0.0_dp
          do I=1,N
             SUM=SUM+DNRM2(N,A(I),LDA)**2
          enddo
          ANORM=SQRT(SUM)
       else if(SPARSE_A.EQ.1 .OR. SPARSE_A.EQ.4) then  ! Sparse matricies: CMISS compressed row
          SUM=0.0_dp
          do I=1,N
             do IJ=ISR_A(I),ISR_A(I+1)-1
                SUM=SUM+A(IJ)**2
             enddo
          enddo
          ANORM=SQRT(SUM)
          
       else
          !          ERROR='>>Sparsity type not implemented'
          !          GOTO 9999
       endif
       
    else
       !        WRITE(ERROR,'(A,I10,A)') '>>Normtype',NORMTYPE,' not supported'
       !        GOTO 9999
    endif
    
  end subroutine CALC_MAT_NORM


!!! ########################################################################## 

      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER :: INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER :: IPIV( * )
      real(dp) :: A( LDA, * ), B( LDB, * )

!     ..
!
!  Purpose
!  =======
!
!  DGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by DGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
      !LOGICAL            LSAME
! HBK COMMENTED
!     EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
! Joyce commented
!      EXTERNAL           DLASWP, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
           LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )  RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A ! X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, &
                   ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, &
                    NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A' ! X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
                   ONE, A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, &
                    A, LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of DGETRS
!
      END SUBROUTINE DGETRS

!!! ########################################################################## 

      double precision function dasum(n,dx,incx)
!
!     takes the sum of the absolute values.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
!
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2)) &
       + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp

      return
      end function dasum


!!! ########################################################################## 
      subroutine  dcopy(n,dx,incx,dy,incy)
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue

      return
      end SUBROUTINE  dcopy

!!! ########################################################################## 
      DOUBLE PRECISION FUNCTION DNRM2 ( N, X, INCX )
!     .. Scalar Arguments ..
      INTEGER                           INCX, N
!     .. Array Arguments ..
      DOUBLE PRECISION                  X( * )
!     ..
!
!  DNRM2 returns the euclidean norm of a vector via the function
!  name, so that
!
!     DNRM2 := sqrt( x'*x )
!
!
!
!  -- This version written on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to DLASSQ.
!     Sven Hammarling, Nag Ltd.
!
!
!     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     .. Local Scalars ..
      INTEGER               IX
      DOUBLE PRECISION      ABSXI, NORM, SCALE, SSQ
!     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
!     ..
!     .. Executable Statements ..
      IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
!
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SSQ   = ONE   + SSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SSQ   = SSQ   +     ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      END IF
!
      DNRM2 = NORM
      RETURN
!
!     End of DNRM2.
!
      END FUNCTION DNRM2
!!! ########################################################################## 
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
      !EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
!     ..
!     .. External Functions ..
     ! INTEGER            ILAENV
      !EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )  RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'DGETRF', M, N, -1)
!      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,IPIV, 1 )
!
!              Compute block row of U.
!
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
                N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                       N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, &
                       A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of DGETRF
!
      END SUBROUTINE DGETRF
!!! ##########################################################################
      LOGICAL          FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME ) &
        RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR. &
            INTA.GE.145 .AND. INTA.LE.153 .OR. &
            INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR. &
            INTB.GE.145 .AND. INTB.LE.153 .OR. &
            INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
      END FUNCTION LSAME
!!! ##########################################################################
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
                        BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
!      LOGICAL            LSAME

!HARI COMMENTED
!     EXTERNAL           LSAME
!     .. External Subroutines ..
!     EXTERNAL           XERBLA
!     .. Intrinsic Functions ..

      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND. &
              ( .NOT.LSAME( TRANSA, 'C' ) ).AND. &
              ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND. &
              ( .NOT.LSAME( TRANSB, 'C' ) ).AND. &
              ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
         ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
        RETURN
!
!     And if  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B' + beta*C
!
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEMM .
!
      END SUBROUTINE DGEMM

!!! ##########################################################################

      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M!abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!  Further Details
!  ===============
!
!  Modified by
!   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
!
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
!
      RETURN
!
!     End of DLASWP
!
      END SUBROUTINE DLASWP

!!! ##########################################################################

      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, JP
!     ..
!     .. External Functions ..
      !INTEGER            IDAMAX
      !EXTERNAL           IDAMAX
!     ..
!     .. External Subroutines ..
      !EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) &
        RETURN
!
      DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF( JP.NE.J ) &
              CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            IF( J.LT.M ) &
              CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         ELSE IF( INFO.EQ.0 ) THEN
!
            INFO = J
         END IF
!
         IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA, &
                      A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
!
!     End of DGETF2
!
      END SUBROUTINE DGETF2

!!! ##########################################################################

      INTEGER FUNCTION ILAENV( ISPEC, NAME, N1, N2, N4 )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME !, OPTS
      INTEGER            ISPEC, N1, N2, N4 !, N3
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!          = 9: maximum size of the subproblems at the bottom of the
!               computation tree in the divide-and-conquer algorithm
!               (used by xGELSD and xGESDD)
!          =10: ieee NaN arithmetic can be trusted not to trap
!          =11: infinity arithmetic can be trusted not to trap
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. External Functions ..
      INTEGER            IEEECK
      !EXTERNAL           IEEECK
!     ..
!     .. Executable Statements ..
!

      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000, &
             1100 ) ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
  100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) &
                 SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
            ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
            ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                  ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                  ( IC.GE.162 .AND. IC.LE.169 ) ) &
                 SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 ) &
                 SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) ) &
        RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
!
      GO TO ( 110, 200, 300 ) ISPEC
!
  110 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
                 C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
               C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
               C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
               C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
               C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
               C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
               C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
               C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
               C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
!
  200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
               C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
               C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
               C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
               C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
               C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
               C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
               C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
               C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
  300 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
               C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
               C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
               C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
               C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
!
  400 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
  500 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  600 CONTINUE 
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
  900 CONTINUE
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
      ILAENV = 25
      RETURN
!
 1000 CONTINUE
!
!     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
!
!!     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 ) 
      END IF
      RETURN
!
 1100 CONTINUE
!
!     ISPEC = 11: infinity arithmetic can be trusted not to trap
!
!!     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 ) 
      END IF
      RETURN
!
!     End of ILAENV
!
      END FUNCTION ILAENV

!!! ##########################################################################

      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
           B, LDB )
!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      !LOGICAL            LSAME
! HBK COMMENED
!     EXTERNAL           LSAME
!     .. External Subroutines ..
      !EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND. &
             ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND. &
              ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND. &
              ( .NOT.LSAME( TRANSA, 'T' ) ).AND. &
              ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND. &
              ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 ) &
        RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT ) &
                          B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT ) &
                          B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B.
!
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT ) &
                       TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT ) &
                       TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' ).
!
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of DTRSM .
!
      END SUBROUTINE DTRSM
!!! ##########################################################################
      SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (preliminary version) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
!
      WRITE( *, FMT = 9999 )SRNAME, INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
           'an illegal value' )
!
!     End of XERBLA
!
      END SUBROUTINE XERBLA
!!! ##########################################################################
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1998
!
!     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  INTEGER
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, &
                        NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
      IEEECK = 1
!
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      IF( ISPEC.EQ.0 ) &
        RETURN
!
      NAN1 = POSINF + NEGINF
!
      NAN2 = POSINF / NEGINF
!
      NAN3 = POSINF / POSINF
!
      NAN4 = POSINF*ZERO
!
      NAN5 = NEGINF*NEGZRO
!
      NAN6 = NAN5*0.0
!
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      RETURN
      END FUNCTION IEEECK
!!! ##########################################################################
!!! ##########################################################################
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
      !EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) ) &
        RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of DGER  .
!
      END SUBROUTINE dger

!####################################################################
      subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end subroutine dscal

!#######################################################
      subroutine  dswap (n,dx,incx,dy,incy)
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue

      return
      end subroutine dswap

!####################################################################

      integer function idamax(n,dx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision :: dx(*),dmax
      integer :: i,incx,ix,n
!
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue

      return
      end function idamax

!##########################################################


end module solve
