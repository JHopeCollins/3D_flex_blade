!subroutine to construct influence matrix and rhs of impermeability
!condition

subroutine matrix(xc,xv,xn,a,b,rhs)

use constants
use calls
use options

!variables
implicit none

!inputs
real, intent(in), dimension(N,N,3) :: xc        !collocation points
real, intent(in), dimension(N,N,3) :: xv        !vortex ring corners
real, intent(in), dimension(N,N,3) :: xn        !panel normals

!outputs
real, intent(out), dimension(N,N) :: a          !influence matrix
real, intent(out), dimension(N,N) :: b          !induced downwash influence matrix
real, intent(out), dimension(N)   :: rhs        !rhs of linear system

!internal variables
integer :: i,j, p,q
integer :: countC,countV
real, dimension(3) :: u,v       !velocities
real, dimension(2,3) :: ut

!initialise
u   = 0.0
v   = 0.0
a   = 0.0
b   = 0.0
rhs = 0.0

!influence matrix arranged in blocks of the blade section at each spanwise location
!scan over collocation points
do 1 j = 1,Ns
do 1 i = 1,Nc
        
        countC = indx(i,j,Nc)
        
        rhs(countC) = -dot_product(uinf , xn(i,j,:))
        
        write(160,*) 'c', countC, 'rhs', rhs(countC)
        
        do 1 q = 1,Ns

                !scan over bound vortex rings
                do 2 p = 1,Nc
                        
                        countV = indx(p,q,Nc)
                        
                        write(150,*) 'c', xc(i,j,1:2), 'v', xv(p,q,1:2), countC,countV
                        
                        ut = vring_trail(xc(i,j,:), xv(p  ,q  ,:), xv(p  ,q+1,:), &
                                                    xv(p+1,q+1,:), xv(p+1,q  ,:))
                        u  = ut(1,:)
                        v  = ut(2,:)
                        
                        a(countC,countV) = a(countC,countV) + dot_product(u , xn(i,j,:))
                        b(countC,countV) = b(countC,countV) + dot_product(v , xn(i,j,:))
                        
                2 continue
                
                countV = indx(Nc,q,Nc)

                write(150,*) 'c', xc(i,j,1:2), countC,countV
                write(150,*) ''
                
                !scan over shed vortex rings
                do 1 p = Nc+1,Nw-1
                        
                        ut = vring_trail(xc(i,j,:), xv(p  ,q  ,:), xv(p  ,q+1,:), &
                                                    xv(p+1,q+1,:), xv(p+1,q  ,:))
                        u  = ut(1,:)
                        v  = ut(2,:)
                        
                        a(countC,countV) = a(countC,countV) + dot_product(u , xn(i,j,:))
                        b(countC,countV) = b(countC,countV) + dot_product(v , xn(i,j,:))
        
1 continue


end subroutine matrix

!--------------------------------------------------------------------------


!influence matrix and rhs of flexible impermeability
!condition

subroutine flex_matrix(xc, xn, xv, dxc, dxn, da, drhs)

use constants
use calls
use options

!variables
implicit none

!inputs
real,    intent(in), dimension(N,N,3) :: xc        !collocation points
real,    intent(in), dimension(N,N,3) :: xv        !vortex ring corners
real,    intent(in), dimension(N,N,3) :: xn        !panel normals
complex, intent(in), dimension(N,N,3) :: dxc       !collocation points
complex, intent(in), dimension(N,N,3) :: dxn       !panel normals

!outputs
complex, intent(out), dimension(N,N) :: da          !influence matrix
complex, intent(out), dimension(N)   :: drhs        !rhs of linear system

!internal variables
integer :: i,j, p,q, r
integer :: countC,countV
real, dimension(3) :: u
complex :: v
complex :: phase0, phase, shift

!initialise
u    = 0.0
v    = 0.0
da   = 0.0
drhs = 0.0

phase0 = -imag*omega*dt
!phase0 = 1.0-exp(-imag*omega*dt)
shift  = exp(imag*omega*dt)

!influence matrix arranged in blocks of the blade section at each spanwise location
!scan over collocation points
do 1 j = 1,Ns
do 1 i = 1,Nc
        
        countC = indx(i,j,Nc)
        
        !u.d(xn), and (dxc/dt).xn
        drhs(countC) = -imag*omega*dot_product(dxc(i,j,:) ,  xn(i,j,:)) &
                                 - dot_product(uinf       , dxn(i,j,:))
        
        do 1 q = 1,Ns

                !scan over bound vortex rings
                do 2 p = 1,Nc
                        
                        countV = indx(p,q,Nc)
                        
                        u = vring(xc(i,j,:), xv(p  ,q  ,:), xv(p  ,q+1,:), &
                                             xv(p+1,q+1,:), xv(p+1,q  ,:))
                        
                        da(countC,countV) = da(countC,countV) + dot_product(u , xn(i,j,:))
                        
                2 continue
                
                !scan over shed vortex rings
                do 1 p = Nc+1,Nw-1
                        
                        phase = phase0*exp(imag*omega* &
                                       0.25*(xv(p  ,q  ,1) + xv(p  ,q+1,1) &
                                           + xv(p+1,q+1,1) + xv(p+1,q  ,1)))
                        
                        u = vring(xc(i,j,:), xv(p  ,q  ,:), xv(p  ,q+1,:), &
                                             xv(p+1,q+1,:), xv(p+1,q  ,:))
                        
                        v = phase*dot_product(u , xn(i,j,:))
                        
                        do 1 r = 1,Nc
                                
                                countV = indx(r,q,Nc)
                                
                                da(countC,countV) = da(countC,countV) + v
        
1 continue


end subroutine flex_matrix
