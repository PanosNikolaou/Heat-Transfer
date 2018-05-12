module tools

  contains

  subroutine printmtx(a)
  use types

    real(kind=rk), dimension (:,:), intent(in) :: a
    integer :: i

        do i = 1, ubound (a,1)
            print *, "row numbers", i, "/", a(i,:)
        end do
    end subroutine printmtx

  subroutine gauss(a,b,x)
  use types

    implicit none

    real (kind=rk), dimension (:,:), intent (in),allocatable :: a
    real (kind=rk), dimension (:), intent (out),allocatable :: x
    real (kind=rk), dimension (:), intent (in),allocatable :: b

    real (kind=rk), dimension (:,:), allocatable :: acopy

    integer :: i, j, n

    if (.not. allocated(a) .or. .not.(allocated(b))) then
        print *, "inputs to gem are not allocated..."
        return
    end if

    if ((ubound(a,1) /= ubound(a,2))) then
        print*, "matrix a is not square..."
        return
    end if

    if ((ubound(a,1) /= ubound(b,1))) then
        print*, "vector b and matrix a have incompatible lengths..."
        return
    end if

    n = ubound(a,1)

    if (.not. allocated (x)) then
        allocate(x(n))
    end if

    allocate(acopy(ubound(a,1), ubound(a,1)+1))
    acopy(:, 1:ubound(a,1)) = a
    acopy(:, ubound(a,1)+1) = b

    do j=2, ubound(a,1)
       do i = j, ubound(a,1)
         acopy(i, :) = acopy(i,:) - acopy(j-1,:)/acopy(j-1,j-1)*acopy(i,j-1)
       end do
    end do

    forall (i=n:1:-1) x(i) = ( acopy(i,n+1) - sum(acopy(i,i+1:n)*x(i+1:n)) ) / acopy(i,i)

    !call printmtx(acopy)

  end subroutine gauss

  subroutine swap_reals(a, b)
  use types
    real (kind=rk), intent(inout) :: a, b
    real (kind=rk) :: tmp

    tmp = a
    a = b
    b = tmp
  end subroutine swap_reals

  subroutine gen_hilbert(n, A)
  use types
    integer, intent(in) :: n
    integer :: i,j

    real (kind=rk), dimension (:,:), intent (out), allocatable :: A

    allocate(A(n,n))
    do i=1,n
        do j=1,n
            A(i,j) = 1.0_rk/(i+j-1)
        end do
    end do
  end subroutine gen_hilbert

  subroutine inverse(a,c,n)
  !============================================================
  ! Inverse matrix
  ! Method: Based on Doolittle LU factorization for Ax=b
  ! Alex G. December 2009
  !-----------------------------------------------------------
  ! input ...
  ! a(n,n) - array of coefficients for matrix A
  ! n      - dimension
  ! output ...
  ! c(n,n) - inverse matrix of A
  ! comments ...
  ! the original matrix a(n,n) will be destroyed
  ! during the calculation
  !===========================================================
  use types
    implicit none
    integer n
    real (kind=rk) ::  a(n,n), c(n,n)
    real (kind=rk) ::  L(n,n), U(n,n), b(n), d(n), x(n)
    real (kind=rk) ::  coeff
    integer i, j, k

    ! step 0: initialization for matrices L and U and b
    L=0.0
    U=0.0
    b=0.0

    ! step 1: forward elimination
    do k=1, n-1
        do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
                do j=k+1,n
                    a(i,j) = a(i,j)-coeff*a(k,j)
                end do
        end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
        L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
        do i=1,j
            U(i,j) = a(i,j)
        end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
        b(k)=1.0
        d(1) = b(1)

        ! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
            d(i)=b(i)
                do j=1,i-1
                    d(i) = d(i) - L(i,j)*d(j)
                end do
        end do

        ! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
            x(i) = d(i)
                do j=n,i+1,-1
                    x(i)=x(i)-U(i,j)*x(j)
                end do
            x(i) = x(i)/u(i,i)
        end do

        ! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
            c(i,k) = x(i)
        end do
        b(k)=0.0
    end do
  end subroutine inverse

! ---------------------- MODULE fseidel.f90 ---------------------------
!Gauss Seidel Method with relaxation

subroutine seidel(crit,n,mat,b,omega,x,residu,iter,rc)
parameter(ITERMAX=500)            ! Maximal number of iterations
parameter(ONE=1.d0,TWO=2.d0,ZERO=0.d0)
  integer crit, n, iter, rc
  REAL*8 mat(n,n),b(n),omega
  REAL*8 x(n),residu(n)
!*====================================================================*
!*                                                                    *
!*  seidel solves the linear system  mat * x = b  iteratively.        *
!*  Here  mat  is a nonsingular  n x n  matrix, b is the right hand   *
!*  side for the linear system and x is the solution.                 *
!*                                                                    *
!*  seidel uses the Gauss Seidel Method with relaxation for a given   *
!*  relaxation coefficient 0 < omega < 2.                             *
!*  If  omega = 1, the standard Gauss Seidel method (without          *
!*  relaxation) is performed.                                         *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Applications:                                                    *
!*   =============                                                    *
!*      Solve linear systems with nonsingular system matrices that    *
!*      satisfy one of the following criteria: row sum criterion,     *
!*      column sum criterion or the criterion of Schmidt and v. Mises.*
!*      Only if at least one of these criteria is satisfied for mat,  *
!*      convergence of the scheme is guaranteed [See BIBLI 11].       *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Input parameters:                                                *
!*   ================                                                 *
!*      crit     integer crit                                         *
!*               select criterion                                     *
!*               =1 : row sum criterion                               *
!*               =2 : column sum criterion                            *
!*               =3 : criterion of Schmidt-v.Mises                    *
!*               other : no check                                     *
!*      n        integer n ( n > 0 )                                  *
!*               size of mat, b and x                                 *
!*      mat      REAL*8   mat(n,n)                                    *
!*               Matrix of the liear system                           *
!*      b        REAL*8 b(n)                                          *
!*               Right hand side                                      *
!*      omega    REAL*8 omega; (0 < omega < 2)                        *
!*               Relaxation coefficient.                              *
!*      x        REAL*8  x(n)                                         *
!*               Starting vector for iteration                        *
!*                                                                    *
!*   Output parameters:                                               *
!*   ==================                                               *
!*      x        REAL*8  x(n)                                         *
!*               solution vector                                      *
!*      residu   REAL*8   residu(n)                                   *
!*               residual vector  b - mat * x; close to zero vector   *
!*      iter     integer iter                                         *
!*               Number of iterations performed                       *
!*      rc       integer return code                                  *
!*               =  0     solution has been found                     *
!*               =  1     n < 1  or omega <= 0 or omega >= 2          *
!*               =  2     improper mat or b or x (not used here)      *
!*               =  3     one diagonal element of mat vanishes        *
!*               =  4     Iteration number exceeded                   *
!*               = 11     column sum criterion violated               *
!*               = 12     row sum criterion violated                  *
!*               = 13     Schmidt-v.Mises criterion violated          *
!*                                                                    *
!*====================================================================*
  REAL*8 tmp, eps;

  iter = 0                        !Initialize iteration counter
  rc = 0

  if (n<1.or.omega<=ZERO.or.omega>=TWO) then
    rc=1
    return
  end if

  eps = 1.d-10

  do i=1, n                       !transform mat so that all
                                          !diagonals equal 1
    if (mat(i,i) == ZERO) then
      rc=3
      return
    end if
    tmp = ONE / mat(i,i)
    do j=1, n
      mat(i,j)= mat(i,j)*tmp
    end do
    b(i) = b(i)*tmp               !adjust right hand side b

  end do


  !check convergence criteria
  if (crit==1) then
     do i = 1, n                  !row sum criterion
       tmp=ZERO
       do j=1,n
         tmp = tmp + dabs(mat(i,j))
       end do
       if (tmp >= TWO) then
         rc=11
         return
       end if
     end do
  else if (crit==2) then
     do j=1, n                    !column sum criterion
       tmp=ZERO
       do i=1,n
         tmp = tmp + dabs(mat(i,j))
       end do
       if (tmp >= TWO) then
         rc=12
     return
       end if
     end do
  else if (crit==3) then
     tmp=ZERO
     do i=1, n
       do j=1, n                  !criterion of Schmidt
         tmp = tmp + mat(i,j)**2  !von Mises
       end do
     end do
     tmp = DSQRT(tmp - ONE)
     if (tmp >= ONE) then
       rc=13
       return
     end if
  end if

  do i=1, n
    residu(i) = x(i)              !store x in residu
  end do

  do while (iter <= ITERMAX)      !Begin iteration

    iter=iter+1

    do i=1, n
      tmp=b(i)
      do j=1, n
        tmp =  tmp - mat(i,j) * residu(j)
      end do
      residu(i) = residu(i) + omega * tmp
    end do

    do i=1, n                     !check break-off criterion
      tmp = x(i) - residu(i)
      if (DABS (tmp) <= eps) then
        x(i) = residu(i)          !If rc = 0 at end of loop
        rc = 0                    !  -> stop iteration
      else
        do j=1, n
          x(j) = residu(j)
        end do
        rc = 4
        goto 10
      end if
    end do
    if (rc == 0) goto 20          !solution found
10 end do                         !End iteration

20 do i=1, n                      !find residual vector
     tmp=b(i)
     do j=1, n
       tmp = tmp - mat(i,j) * x(j)
     end do
     residu(i) = tmp
   end do

  return

end

      subroutine solve_tridiag(a,b,c,d,x,n)
      implicit none
!    a - sub-diagonal (means it is the diagonal below the main diagonal)
!    b - the main diagonal
!    c - sup-diagonal (means it is the diagonal above the main diagonal)
!    d - right part
!    x - the answer
!    n - number of equations

        integer,intent(in) :: n
        real(8),dimension(n),intent(in) :: a,b,c,d
        real(8),dimension(n),intent(out) :: x
        real(8),dimension(n) :: cp,dp
        real(8) :: m
        integer i

! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
         enddo
! initialize x
         x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do

    end subroutine solve_tridiag

! ----------------------- END fseidel.f90 ----------------------------

end module tools
