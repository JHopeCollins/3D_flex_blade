!load calculation for 3d thin wing under steady free stream

subroutine loads(x, g, lift, drag)

use constants
use calls
use options

!variables
implicit none

!inputs
real, intent(in), dimension(N,N,3) :: x
real, dimension(N,N)   :: g

!outputs
real, intent(out) :: lift
real, intent(out) :: drag

!internal variables
real, dimension(N,N) :: dl, dd
real :: dy
integer :: i, j

!initialise
lift = 0.0
drag = 0.0
dl = 0.0
dd = 0.0
dy = 0.0

do 1 j = 1,Ns
        
        dy =      x(1,j+1,2) - x(1,j,2)
        dy = dy + x(2,j+1,2) - x(2,j,2)
        dy = dy /2.0
        
        dl(1,j) = g(1,j)*dy

1 continue

do 2 j = 1,Ns
do 2 i = 2,Nc
        
        dy =      x(i  ,j+1,2) - x(i  ,j,2)
        dy = dy + x(i+1,j+1,2) - x(i+1,j,2)
        dy = dy /2.0
        
        dl(i,j) = (g(i,j)- g(i-1,j))*dy

2 continue

lift = sum(dl(1:Nc,1:Ns))
drag = sum(dd(1:Nc,1:Ns))


end subroutine loads
