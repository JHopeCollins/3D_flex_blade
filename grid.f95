!construct grid point (panel corners) array, and panel geometry
subroutine grid(x,xc,xv,xn)

use constants
use calls
use options

!variables
implicit none

!outputs
real, intent(out), dimension(N,N,3) :: x        !panel corners
real, intent(out), dimension(N,N,3) :: xc       !collocation points
real, intent(out), dimension(N,N,3) :: xv       !vortex ring corners
real, intent(out), dimension(N,N,3) :: xn       !panel normal vectors

!internal variables
integer :: i,j
real :: xspan   !local spanwise coordinate
real :: xchord  !local chordwise coordinate
real :: chord   !local chord

real, dimension(3) :: Ak, Bk    !panel diagonals

!initialise
x  = 0.0
xc = 0.0
xv = 0.0
xn = 0.0
Ak = 0.0
Bk = 0.0
xchord = 0.0
xspan  = 0.0
chord  = 0.0

!panel corners
x(1,:,2) = linspace(Ns+1, -clip*AR/2.0, clip*AR/2.0)
!x(1,:,2) = cosspace(Ns+1, -clip*AR/2.0, clip*AR/2.0)

do 1 j = 1,Ns+1
        
        x(2:Nc+1,j,2)    = x(1,j,2)
        
        !x(1:Nc+1:Nc,j,1) = flatLEblade(x(1,j,2)*2.0/AR)
        !x(1:Nc+1:Nc,j,1) = flatTEblade(x(1,j,2)*2.0/AR)
        !x(1:Nc+1:Nc,j,1) = sweptblade( x(1,j,2)*2.0/AR)
        x(1:Nc+1:Nc,j,1) = ellipsoid(  x(1,j,2)*2.0/AR)
        
        x(:,j,1) = linspace(Nc+1, x(1,j,1), x(Nc+1,j,1))
        !x(:,j,1) = cosspace(Nc+1, x(1,j,1), x(Nc+1,j,1))
        
        do 1 i = 1,Nc+1
                
                !x(i,j,3) = x(Nc+1,j,1) - x(1,j,1)
                !x(i,j,3) = x(i   ,j,1) - x(1,j,1) / x(i,j,3)
                !x(i,j,3) = NACA2c(h,m,x(i,j,3))
                !x(i,j,3) = (x(Nc+1,j,1) - x(1,j,1))*x(i,j,3)

1 continue


!panel normals
do 2 j = 1,Ns
do 2 i = 1,Nc
        
        !panel diagonals
        Ak(1) = x(i+1,j+1,1) - x(i,j,1)
        Ak(2) = x(i+1,j+1,2) - x(i,j,2)
        Ak(3) = x(i+1,j+1,3) - x(i,j,3)
        
        Bk(1) = x(i,j+1,1) - x(i+1,j,1)
        Bk(2) = x(i,j+1,2) - x(i+1,j,2)
        Bk(3) = x(i,j+1,3) - x(i+1,j,3)
        
        !cross product diags
        xn(i,j,1) = Ak(2)*Bk(3) - Ak(3)*Bk(2)
        xn(i,j,2) = Ak(3)*Bk(1) - Ak(1)*Bk(3)
        xn(i,j,3) = Ak(1)*Bk(2) - Ak(2)*Bk(1)
        
        xn(i,j,:) = xn(i,j,:)/norm2(xn(i,j,:))

2 continue


!collocation points at panel centre span, 3/4 chord
do 3 j = 1,Ns
do 3 i = 1,Nc
        
        xc(i,j,1) = 0.125*( x(i  ,j,1) + x(i  ,j+1,1) ) + &
                    0.375*( x(i+1,j,1) + x(i+1,j+1,1) )
        
        xc(i,j,2) = 0.25 *( x(i  ,j,2) + x(i+1,j  ,2)   + &
                            x(i,j+1,2) + x(i+1,j+1,2) )
        
        xc(i,j,3) = 0.125*( x(i  ,j,3) + x(i  ,j+1,3) ) + &
                    0.375*( x(i+1,j,3) + x(i+1,j+1,3) )

3 continue


!bound vortex ring corners. leading vortex line along panel 1/4 chord
do 4 j = 1,Ns+1
do 10 i = 1,Nc
        
        xv(i,j,1) = 0.75*x(i,j,1) + 0.25*x(i+1,j,1)
        xv(i,j,2) = 0.75*x(i,j,2) + 0.25*x(i+1,j,2)
        xv(i,j,3) = 0.75*x(i,j,3) + 0.25*x(i+1,j,3)
        
        write (141,*) xv(i,j,:)

10 continue
        write (141,*) ''
4 continue


!shed vortex ring corners. vortex lines convected with mean flow
do 5 j = 1,Ns+1
do 11 i = Nc+1,Nw
        
        xv(i,j,1) = xv(i-1,j,1) + dt*uinf(1)
        xv(i,j,2) = xv(i-1,j,2) + dt*uinf(2)
        xv(i,j,3) = xv(i-1,j,3) + dt*uinf(3)
        
        write (141,*) xv(i,j,:)

11 continue
        write (141,*) ''
5 continue

do 99 j = 1,Ns+1
do 99 i = 1,Nc+1
        write(142,*) x(i,j,:)
99 continue



end subroutine grid


!---------------------------------------------------------------------------


!change grid points (panel corners) due to blade flexing
subroutine flex_grid(x, dx, dxc, dxv, dxn)

use constants
use calls
use options

!variables
implicit none

!outputs
real,    intent(in),  dimension(N,N,3) :: x         !panel corners
complex, intent(out), dimension(N,N,3) :: dx        !panel corners
complex, intent(out), dimension(N,N,3) :: dxc       !collocation points
complex, intent(out), dimension(N,N,3) :: dxv       !vortex ring corners
complex, intent(out), dimension(N,N,3) :: dxn       !panel normal vectors

!internal variables
integer :: i,j
real :: phase
real, dimension(3) ::  Ak,  Bk    !panel diagonals
complex, dimension(3) :: dAk, dBk    !panel diagonals
real, dimension(2) :: no,ta

!initialise
dx  = 0.0
dxc = 0.0
dxv = 0.0
dxn = 0.0
Ak  = 0.0
Bk  = 0.0
dAk = 0.0
dBk = 0.0


do 1 j = 1,Ns+1
        
        do 1 i = 1,Nc+1
                
                phase = omega*x(i,j,1)
                
                call NACA2c_vec(h,m,x(i,j,1),no,ta)
                
                dx(i,j,1) = cmplx(1.0,0.0)*no(1)
                dx(i,j,2) = 0.0
                dx(i,j,3) = cmplx(1.0,0.0)*no(2)
                dx(i,j,:) = dx(i,j,:) * exp(imag*phase)

1 continue


!change in panel normals
! d(a x b) = (da x b) + (a x db)
do 2 j = 1,Ns
do 2 i = 1,Nc
        
        !panel diagonals
        Ak(1) = x(i+1,j+1,1) - x(i,j,1)
        Ak(2) = x(i+1,j+1,2) - x(i,j,2)
        Ak(3) = x(i+1,j+1,3) - x(i,j,3)
        
        Bk(1) = x(i,j+1,1) - x(i+1,j,1)
        Bk(2) = x(i,j+1,2) - x(i+1,j,2)
        Bk(3) = x(i,j+1,3) - x(i+1,j,3)
        
        !panel diagonals
        dAk(1) = dx(i+1,j+1,1) - dx(i,j,1)
        dAk(2) = dx(i+1,j+1,2) - dx(i,j,2)
        dAk(3) = dx(i+1,j+1,3) - dx(i,j,3)
        
        dBk(1) = dx(i,j+1,1) - dx(i+1,j,1)
        dBk(2) = dx(i,j+1,2) - dx(i+1,j,2)
        dBk(3) = dx(i,j+1,3) - dx(i+1,j,3)
        
        !cross product diags
        dxn(i,j,1) = dAk(2)*Bk(3)-dAk(3)*Bk(2) + Ak(2)*dBk(3)-Ak(3)*dBk(2)
        dxn(i,j,2) = dAk(3)*Bk(1)-dAk(1)*Bk(3) + Ak(3)*dBk(1)-Ak(1)*dBk(3)
        dxn(i,j,3) = dAk(1)*Bk(2)-dAk(2)*Bk(1) + Ak(1)*dBk(2)-Ak(2)*dBk(1)
        
        dxn(i,j,:) = dxn(i,j,:)/norm2c(3,dxn(i,j,:))

2 continue


!change in collocation points at panel centre span, 3/4 chord
do 3 j = 1,Ns
do 3 i = 1,Nc
        
        dxc(i,j,1) = 0.125*( dx(i  ,j,1) + dx(i  ,j+1,1) ) + &
                     0.375*( dx(i+1,j,1) + dx(i+1,j+1,1) )
        
        dxc(i,j,2) = 0.25 *( dx(i  ,j,2) + dx(i+1,j  ,2)   + &
                             dx(i,j+1,2) + dx(i+1,j+1,2) )
        
        dxc(i,j,3) = 0.125*( dx(i  ,j,3) + dx(i  ,j+1,3) ) + &
                     0.375*( dx(i+1,j,3) + dx(i+1,j+1,3) )

3 continue


!change in bound vortex ring corners. leading vortex line along panel 1/4 chord
do 4 j = 1,Ns+1
do 4 i = 1,Nc
        
        dxv(i,j,1) = 0.75*dx(i,j,1) + 0.25*dx(i+1,j,1)
        dxv(i,j,2) = 0.75*dx(i,j,2) + 0.25*dx(i+1,j,2)
        dxv(i,j,3) = 0.75*dx(i,j,3) + 0.25*dx(i+1,j,3)

4 continue


!change in shed vortex ring corners. vortex lines convected with mean flow
do 5 j = 1,Ns+1
do 5 i = Nc+1,Nw
        
        dxv(i,j,1) = 0.0
        dxv(i,j,2) = 0.0
        dxv(i,j,3) = 0.0

5 continue


end subroutine flex_grid
