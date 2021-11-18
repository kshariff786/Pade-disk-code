!------------------------------------------------------------------------------81

module physical_constants
   implicit none
   real(8), parameter :: Gconst = 6.673d-8 ! cm^3 g^-1 s^-2

   ! Stefan-Boltzmann:
   real(8), parameter :: sigma_SB = 5.6705d-5 ! erg s^-1 cm^-2 K^-4
   ! Boltzmann:
   real(8), parameter :: kB       = 1.3807d-16 ! erg K^-1
   real(8), parameter :: Msolar   = 1.989d33 ! gm
   real(8), parameter :: Rsolar   = 69.57d9  ! cm
   real(8), parameter :: AU       = 1.496d13 ! cm
   real(8), parameter :: parsec   = 3.0857d18 ! cm
   real(8), parameter :: kilometer_per_sec = 1000.d0*100.d0 ! cm/sec
   real(8), parameter :: Rgas     = 3.5871d7  ! cm^2 s^-2 K^-1
   real(8), parameter :: year     = 3.16d7  ! s
   ! Mass of atomic hydrogen:
   real(8), parameter :: mH       = 1.6735d-24 ! gm
   real(8), parameter :: mH2      = 2.d0 * mH  ! gm
end module physical_constants

!------------------------------------------------------------------------------81

use physical_constants
implicit none
real(8) :: l, tau_eddy, pi
pi = 4.0d0 * ATAN(1.0d0)

write (6, 1)
1 format (' enter scale in AU--->', $)
read (5, *) l

print *, ' tau_eddy = ', tau_eddy(l*AU, pi)
stop
end

!------------------------------------------------------------------------------81

! Eddy turnover time at scale l
real(8) function tau_eddy(l, pi)

use physical_constants
implicit none
real(8) :: l, pi

real(8), parameter :: pi = 
! Velocity dispersion:
real(8), parameter :: sigma_1D = 0.13d0 * kilometer_per_sec
! Size of the dense core:
real(8), parameter :: lambda_core = 0.1d0 * parsec
! Kolmogorov constant for 1D spectra:
real(8), parameter :: C1 = 0.5d0

! Rate of dissipation:
real(8) :: epsilon

epsilon = (2.d0/3.d0*sigma_1D**2 / C1)**1.5d0 * 2.d0 * pi / lambda_core
tau_eddy = (l**2 / epsilon)**(1.d0/3.d0)
return
end function tau_eddy

!------------------------------------------------------------------------------81
