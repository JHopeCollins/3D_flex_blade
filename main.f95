!panel method model for 3d thin flexible blade with small deformation and NACA 2 digit camberline and
!prescribed wake

program main

use constants
use calls
use options

!variables
implicit none

real, dimension(N,N,3) :: x             !grid points (panel corners)
real, dimension(N,N,3) :: xc            !panel collocation points
real, dimension(N,N,3) :: xv            !vortex ring corners
real, dimension(N,N,3) :: xn            !panel normals

complex, dimension(N,N,3) :: dx         !grid points (panel corners)
complex, dimension(N,N,3) :: dxc        !panel collocation points
complex, dimension(N,N,3) :: dxv        !vortex ring corners
complex, dimension(N,N,3) :: dxn        !panel normals

real, dimension(N,N) :: a!, adum         !influence matrix
real, dimension(N,N) :: b               !induced drag influence matrix
real, dimension(N)   :: rhs             !free stream influence

complex, dimension(N,N) :: da, dadum    !influence matrix
complex, dimension(N)   :: drhs         !free stream influence

real, dimension(N,N)    :: g            !bound vorticity distribution
complex, dimension(N,N) :: dg           !bound vorticity distribution
real :: lift                            !lift on blade
real :: drag                            !wake induced drag

integer :: i, j, s

!lapack variables
integer :: info1, info2, lda = (Ns*Nc)
integer :: infoeig
integer, dimension(N)   :: ipiv
complex, dimension(5*N) :: CWORK         !cgeev workspace
real, dimension(2*N)    :: RWORK         !cgeev workspace
!real, dimension(5*N)    :: SWORK         !cgeev workspace

!eigen variables
complex, dimension(N)   :: deigval       !eigenvalues of da
complex, dimension(N,N) :: deigvecL      !left  eigenvector i of da is deigvecL(:,i)
complex, dimension(N,N) :: deigvecR      !right eigenvector i of da is deigvecR(:,i)

!real, dimension(N)    :: eigval_real     !real part of eigenvalues of a
!real, dimension(N)    :: eigval_imag     !imag part of eigenvalues of a
!complex, dimension(N) :: eigval          !eigenvalues of a
!real, dimension(N,N)  :: eigvecL         !left  eigenvector i of a is eigvecL(:,i)
!real, dimension(N,N)  :: eigvecR         !right eigenvector i of a is eigvecR(:,i)

!initialise
g = 0.0

!----------------------------------------------------------
!steady

!panel grid
call grid(x, xc, xv, xn)


!influence matrix and rhs
call matrix(xc, xv, xn, a, b, rhs)

!adum = a

!call sgeev('V', 'V', lda, a, N, eigval_real, eigval_imag, eigvecL, N, eigvecR, N, SWORK, 5*N, infoeig)
!print *, 'steady eig info', infoeig
!eigval = cmplx(eigval_real, eigval_imag)

!do 120 i = lda,1, -1
!        print *, i, eigval_real(i), eigval_imag(i)
!120 continue

!do 100 i = 1,lda
 !       eigval(i) = dumeigval(maxloc(dumeigval))




!a = adum

!call sgetrf(     lda, lda, a, N, ipiv,         info1)
!call sgetrs('n', lda, 1,   a, N, ipiv, rhs, N, info2)
!print *, 'steady info1', info1
!print *, 'steady info2', info2
!!g(1:Nc,1:Ns) = rvec2mat(lda, Nc, Ns, rhs(1:lda))
!do 90 j = 1,Ns
!do 90 i = 1,Nc
!        g(i,j) = rhs(indx(i,j,Nc))
!90 continue
!
!open(unit=3, file='Sgam.dat', action='write', status='replace')
!do 91 j = 1,Ns
!        write(3,*) xc(1,j,1:2), (g(1,j)), g(1,j)
!do 92 i = 2,Nc
!        write(3,*) xc(i,j,1:2), g(i,j), g(i,j)-g(i-1,j)
!92 continue
!91 continue
!
!
!open(unit=4, file='Slefteigvec.dat', action='write', status='replace')
!do 96 s = 1,10
!do 93 j = 1,Ns
!do 93 i = 1,Nc
!        write(1000+s,*) xc(i,j,1:2), eigvecL(indx(i,j,Nc),s)
!93 continue
!96 continue


call loads(x, g, lift, drag)


!----------------------------------------------------------
!unsteady


!change in panel grid
call flex_grid(x, dx, dxc, dxv, dxn)


!influence matrix and rhs for flexible oscillations
call flex_matrix(xc, xn, xv, dxc, dxn, da, drhs)

dadum = da

call cgeev('V', 'V', lda, da, N, deigval, deigvecL, N, deigvecR, N, CWORK, 5*N, RWORK, infoeig)
print *, 'flex eig info', infoeig

da = dadum

call cgetrf(     lda, lda, da, N, ipiv,          info1)
call cgetrs('n', lda, 1,   da, N, ipiv, drhs, N, info2)
print *, 'flexible info1', info1
print *, 'flexible info2', info2
do 80 j = 1,Ns
do 80 i = 1,Nc
        dg(i,j) = drhs(indx(i,j,Nc))
80 continue


open(unit=2, file='gam.dat', action='write', status='replace')
do 97 j = 1,Ns
        write(2,*) xc(1,j,1:2), real(dg(1,j)), aimag(dg(1,j)), real(dg(1,j)), aimag(dg(1,j))
do 98 i = 2,Nc
        write(2,*) xc(i,j,1:2), real(dg(i,j)), aimag(dg(i,j)) &
                , real(dg(i,j)-dg(i-1,j)), aimag(dg(i,j)-dg(i-1,j))
98 continue
97 continue


do 99 s = 1,10
do 99 j = 1,Ns
do 99 i = 1,Nc
        write(2000+s,*) xc(i,j,1:2), real(deigvecL(indx(i,j,Nc),s)), aimag(deigvecL(indx(i,j,Nc),s))
99 continue



end program main
