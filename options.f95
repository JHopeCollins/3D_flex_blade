!set parameters for flexible 3d blade

module options

use constants

!variables
implicit none

!options
!blade geometry
real, parameter :: rootchord =  1.0      !chord length at blade root
real, parameter :: tipchord  =  1.0      !chord length at blade tip
real, parameter :: leadaxis  =  0.5      !relative size of streamwise axis of le ellipsoid
real, parameter :: AR        =  4.0      !aspect ratio relative to root chord
real, parameter :: clip      =  0.8      !clipped wing tips for ellipsoid

!camberline geometry
real, parameter :: h = 0.00             !maximum camber height
real, parameter :: m = 0.40             !maximum camber location

!freestream velocity
real, parameter :: alpha = 0.0*deg2rad
real, parameter, dimension(3) :: uinf = [cos(alpha) , 0.0 , sin(alpha)]

!flexing wavenumbers
real, parameter :: lamda = 2.0                  !wavelength relative to chord
real, parameter :: k     = pi/lamda             !chordwise reduced frequency
real, parameter :: omega = 2.0*k                !frequency / wavenumber

!discretisation values
integer, parameter :: Nwp = 10                          !number of wake wavelengths
integer, parameter :: res = 20                          !wake panels/wavelength
real,    parameter :: panelAR = 4.0                     !panel aspect ratio
integer, parameter :: Nc = int(res/lamda)               !number of chordwise panels
integer, parameter :: Ns = Nc*ceiling(AR/panelAR)   !number of spanwise panels
integer, parameter :: Nw = Nc + Nwp*res                 !number of last wake panel
real,    parameter :: dt  = lamda/res

!derived values
real, parameter    :: span = rootchord*AR    !blade spanwise length
real, parameter    :: trailaxis =  1.0-leadaxis

end module  options
