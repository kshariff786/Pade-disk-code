!----------------------------------------------------------------------------------85

subroutine hydrostatic_test

! Hydrostatic in z with only one grid point in r and phi.  There is no azimuthal
! velocity.

!use hydrostatic_parameters
use grid, only: nz, nr, nphi, ndof, rmin, rmax, zmin, zmax, phi_min, phi_max, zgrid, rgrid
use q_array
! In case the we want to run isothermally:
use thermal_parameters, only: ci_squared_initial
use boundary_condition_types
use logical_units
use activate_routines
use basic_state
implicit none

! Locals:
external rhs
real(8), allocatable, dimension(:) :: rho_bar, p_bar

! For domain (read from namelist):
real(8) :: zmax_over_H

! For thermal parameters:
real(8) :: gamma, gm1
logical :: isothermal

! For gravity set-up: 
integer :: gravity_flag, i_thin, i_no_radial
real(8) :: GM

integer :: istep, nsteps, profiles_interval
real(8) :: cfl, t, dt, Omega, r0
integer :: ir, iz, ir_mid

real(8) :: eps, sigma, z0

! Isothermal and adiabatic sound speeds:
real(8) :: ci, ca

real(8) :: rho0, H
logical :: suppress_z_derivatives_when_nz_not_1 = .false., periodic_z = .false., &
     use_supplied_dt = .false.

! For restart if needed:
logical :: restart
integer :: istep_of_restart_file, file_version
real(8) :: t_of_restart_file

integer :: rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced

logical :: apply_pade_filter, filter_relative_to_basic_state = .true.
character(3) :: eps_or_tau
real(8) :: eps_filter, tau_filter

logical :: apply_artificial_pressure
real(8) :: vke, dt_old, C_ap

logical, parameter :: stretched_r = .false., stretched_z = .false.
real(8), parameter :: r0_unif = 0.d0, z0_unif = 0.d0
integer, parameter :: nr_u = 0, nz_u = 0
logical :: negative_density

namelist /hydrostatic_test_input/ nz, cfl, nsteps, profiles_interval, &
     apply_pade_filter, eps_or_tau, eps_filter, tau_filter, eps, sigma, z0, zmax_over_H, &
     apply_artificial_pressure, C_ap, isothermal

open (unit = lun_general_purpose, file = 'input_file', form = 'formatted', &
      status = 'old')
read (lun_general_purpose, nml = hydrostatic_test_input)
close (lun_general_purpose)

allocate(rho_bar(nz), p_bar(nz))

gamma      = 1.4d0
gm1        = gamma - 1.d0

! For units can set three things to unity:
GM   = 1.0d0
H    = 1.0d0
rho0 = 1.0d0

! This should be considered a parameter:
r0      = 1.0d0*H

Omega = SQRT(GM/r0**3)
ci    = H*Omega/SQRT(2.d0)

ca    = SQRT(gamma) * ci   ! adiabatic sound speed

zmin    = -zmax_over_H*H
zmax    =  zmax_over_H*H
rmin    =  r0
rmax    =  r0
phi_min = 0.0d0
phi_max = 0.0d0
nr      = 1
nphi    = 1

gravity_flag = 1
i_thin       = 1 ! thin disk gravity which is linear in z:
! We don't have any azimuthal velocity so don't put any radial gravity.
i_no_radial = 1
restart = .false.
file_version = 2

rmin_BC   = null
rmax_BC   = null
!zmin_BC   = non_reflective
!zmax_BC   = non_reflective
ibalanced = 1
zmin_BC = zero_normal_momentum
zmax_BC = zero_normal_momentum

print *, ' isothermal = ', isothermal
print *, ' ci = ', ci, ' ca = ', ca

! These three set-up calls should be in every application subroutine:
call set_up_domain_mesh_and_partition(rmin, rmax, zmin, zmax, phi_min, phi_max, nr, nz, nphi, &
     suppress_z_derivatives_when_nz_not_1, periodic_z, &
     stretched_r, stretched_z, r0_unif, nr_u, z0_unif, nz_u)     
call set_up_thermal_parameters(gamma, isothermal)
call set_up_boundary_conditions(rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced)

call activate_gravity(gravity_flag, i_thin, i_no_radial, GM, .false., 0.0d0)
if (apply_pade_filter) then
   call activate_pade_filter(apply_pade_filter, eps_or_tau, eps_filter, &
        tau_filter, filter_relative_to_basic_state)
end if

if (apply_artificial_pressure) then
   call activate_artificial_pressure (C_ap)
end if

if (restart) then
   call read_restart(file_version, istep_of_restart_file, t_of_restart_file)
end if

! Needed for the isothermal option.  ci_squared_initial sits in
! module thermal_parameters.
if (isothermal) then
   do iz = 1, nz
      do ir = 1, nr
         ci_squared_initial(ir, iz) = ci**2
      end do
   end do
end if
   
#ifdef debug_print
   print *, ' Finished initializating ci_squared_initial'
#endif

call initial_condition

#ifdef debug_print
   print *, ' Returned from initial condition'
#endif
      
call store_basic_state(filter_relative_to_basic_state)   
call output_gravity_profile_at_mid_radius   

#ifdef debug_print
   print *, ' Finished hydrostatic initial condition'
#endif

t = 0.0d0
call output
open (unit = 101, file = 'vke.dat', form = 'formatted', status = 'unknown')
do istep = 1, nsteps
   call rk4(rhs, q, cfl, t, dt, use_supplied_dt, negative_density)
   if (negative_density) call terminate_with_save(1, istep-1, t)
   
   print *, ' finished istep = ', istep, ' t = ', t, ' dt = ', dt

   call vertical_kinetic_energy (q, vke)
   write (101, "(2(1x, e12.5))") t, vke
   
   if (mod(istep, profiles_interval) .eq. 0) then
      call output
   end if

   dt_old = dt
   if (istep .gt. 1) then
      if (dt .lt. 0.1d0*dt_old) then
         print *, ' dt decreased by a factor of 10 in one step'
         go to 100
      end if
   end if
end do

100 continue
close(101)
call terminate_with_save(0, istep-1, t)

contains

!----------------------------------------------------------------------------------85

subroutine output

use dof_indices
implicit none
character(40), dimension(5) :: filename
integer       :: iz, ifile
real(8)       :: rho, uz, p, ci

write (filename(1), "('rho_t=',   i3.3, f0.4, '.dat')") int(t), t - int(t)
write (filename(2), "('uz_t=',    i3.3, f0.4, '.dat')") int(t), t - int(t)
write (filename(3), "('p_t=',     i3.3, f0.4, '.dat')") int(t), t - int(t)
write (filename(4), "('rhop_t=',  i3.3, f0.4, '.dat')") int(t), t - int(t)

do ifile = 1, 4
   open (unit = 10+ifile, file = filename(ifile), form = 'formatted', status = 'unknown')
end do

do iz = 1, nz
   rho       = q(1,1, iz, irho)
   uz        = q(1,1, iz, zmom) / rho

   if (isothermal) then
      p = rho * ci**2
   else
      p = q(1,1, iz, ener) * gm1
   end if
   
   write (11, "(2(1x, e13.5e3))") zgrid(iz), rho
   write (12, "(2(1x, e13.5e3))") zgrid(iz), uz
   write (13, "(2(1x, e13.5e3))") zgrid(iz), p
   write (14, "(2(1x, e13.5e3))") zgrid(iz), rho - rho_bar(iz)     

   if (q(1, 1, iz, amom) .ne. 0.d0) then
      print *, ' angular momentum .ne. 0 for iz = ', iz
      stop
   end if

   if (q(1, 1, iz, rmom) .ne. 0.d0) then
      print *, ' radial momentum .ne. 0 for iz = ', iz
      stop
   end if
end do
do ifile = 1, 4
   close(10+ifile)
end do
end subroutine output

!~~~~~~~~~~~~~~~
subroutine initial_condition

use dof_indices, only: irho, rmom, amom, zmom, ener
use grid
use gravity
use thermal_parameters
implicit none

! Local:
integer :: ir, iphi, iz
real(8) :: rho_provisional, c2, N2, pi
real(8), dimension(nz) :: dpdz, u_pert, theta, d_theta_dz, dpdz_exact

real(8) :: dpdz_err, rho_err, c2_err

real(8) :: z_shift

! Relative density perturbation in the notation of Liepmann and Roshko:
real(8):: s_tilde

! For quantities after acoustic perturbation has been added:
real(8) :: rho1, p1, uz1

do iz = 1, nz
   rho_provisional = rho0 * EXP(-zgrid(iz)**2 / H**2)
   p_bar(iz)       = rho_provisional * ci**2
   dpdz_exact(iz)  = p_bar(iz) * (-2.d0 * zgrid(iz) / H**2) 
end do

call pade_diff_z(nr*nphi, p_bar, dpdz)

! Set rho to satisfy numerical hydrostatic balance:
do iz = 1, nz
   ! Required since gz = 0 at the midplane:
   if (gz(1, iz) .ne. 0.d0) then
      rho_bar(iz) = dpdz(iz)/gz(1, iz)
    else
      rho_bar(iz) = rho0
   end if
   print *, ' iz = ', iz, ' rho_bar = ', rho_bar(iz), ' dpdz = ', dpdz(iz), ' rho0 = ', rho0
end do
print *, ' enter anything to continue'
read (5, *)

open (unit = 1, file = 'dp_dz_error.dat', form = 'formatted', status = 'unknown')
open (unit = 2, file = 'c2_error.dat',    form = 'formatted', status = 'unknown')
do iz = 1, nz
   dpdz_err = ABS(dpdz(iz) - dpdz_exact(iz))
   rho_err  = 1.d0 / gz(1, iz) * dpdz_err
   c2_err = - gamma * p_bar(iz) * rho_err / rho_bar(iz)**2
   write (1, "(2(1x, e12.5))") zgrid(iz), dpdz_err
   write (2, "(2(1x, e12.5))") zgrid(iz), c2_err
end do
close (1)
close (2)

print *, ' adding perturbation'

do iz = 1, nz
   do iphi = 1, nphi
      do ir = 1, nr
         print *, ' ir = ', ir, ' iphi = ', iphi, ' iz = ', iz
         z_shift = zgrid(iz) - z0
         s_tilde = eps*exp(-z_shift**2/sigma**2)
         if (isothermal) then
            uz1 = ci * s_tilde
         else
            uz1 = ca * s_tilde
         end if
         rho1 = rho_bar(iz) * (1.d0 + s_tilde)
         p1   =   p_bar(iz) * (1.d0 + s_tilde)
         q(ir, iphi, iz, irho) = rho1
         q(ir, iphi, iz, rmom) = 0.d0
         q(ir, iphi, iz, amom) = 0.d0
         q(ir, iphi, iz, zmom) = rho1 * uz1

         ! In case we are adiabatic:
         q(ir, iphi, iz, ener) = p1 / gm1
      end do
   end do
end do

print *, ' Finished adding perturbation'

! Calculate the Brunt-Vaisala frequency:
do iz = 1, nz
   c2 = gamma * p_bar(iz) / rho_bar(iz)
   ! This is proportional to the potential temperature:
   theta(iz) = c2 * p_bar(iz)**(gamma/gm1)
end do

print *,' calling pade_diff_bundle'
call pade_diff_z(1, theta, d_theta_dz)
print *, ' returned from pade_diff_bundle'

open (unit = 1, file = 'Brunt_Vaisala_Squared.dat', form = 'formatted', &
     status = 'unknown')
do iz = 1, nz
   ! Square of the Brunt-Vaisala frequency:
   N2 = gz(1, iz) / theta(iz) * d_theta_dz(iz)
   write (1, "(2(1x, e12.5))") zgrid(iz), N2
end do
close (1)

return
end subroutine initial_condition
!~~~~~~~

end subroutine hydrostatic_test

!----------------------------------------------------------------------------------85

subroutine vertical_kinetic_energy (q, vke)

use partition_data
use grid
use dof_indices
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: vke

! Local:
real(8) :: integrand_1, integrand_2
integer :: iz

vke = 0.0d0
do iz = 1, nz-1
   integrand_1 = q(1, 1, iz,   zmom)**2 / q(1, 1, iz,   irho)*Jz(iz)
   integrand_2 = q(1, 1, iz+1, zmom)**2 / q(1, 1, iz+1, irho)*Jz(iz+1)
   vke = 0.5d0 * (integrand_1 + integrand_2)
end do
vke = 0.5d0*vke

end subroutine vertical_kinetic_energy

!----------------------------------------------------------------------------------85





