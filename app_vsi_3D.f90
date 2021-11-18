!----------------------------------------------------------------------------------85

subroutine vsi_3D

use dof_indices, only: irho, amom, zmom, rmom, ener
use grid
use boundary_condition_types 
use thermal_parameters, only: ci_squared_initial, gm1, d_ci_dr_inner, d_ci_dr_outer
use gravity, only: gz, gr
use q_array, only: q
use logical_units
use control_parameters, only: n_words_allocated, restart_file, save_file
use math_constants, only: pi
use partition_data
use viscosity_types
use basic_state
use artificial_pressure_module, only: iphi_max_ap, p_art, dil ! For outputting the artificial pressure
use cpu_timing_module
use activate_routines
#ifdef mpi_code
   use mpi
#endif
implicit none

external rhs ! name of rhs subroutine passed to rk4.
logical :: restart
integer :: istep_of_restart_file, file_version
real(8) :: t_of_restart_file

integer :: istep, nsteps, istep0
real(8) :: t, cfl, dt, dt_previous, dt_min
real(8) :: H0, rho0, c0, r0
real(8) :: zmax_over_H0, H0_over_r0, rsize_over_H0, lambda_r, V_Kep0, &
     Omega0, T_orbital0, phi_max_over_pi

! Pade filter stuff:
logical :: apply_pade_filter, filter_relative_to_basic_state =.true.
character(3) :: eps_or_tau
real(8) :: eps_filter, tau_filter

! For thermal set-up:
real(8) :: gamma      = 1.4d0
logical :: isothermal

! For gravity set-up
real(8) :: GM
integer :: gravity_flag, i_thin, i_no_radial

real(8) :: pexp, qexp, c_exponent, H_exponent ! exponents for rho, temp, etc.
! Formula temps.:
real(8) :: Omega, rho_temporary, term1, term2, term3, R_spherical, arg

integer :: iz, ir, iphi, ier
logical :: output_profiles, output_profiles_now
integer :: tecplot_interval, profiles_interval, fluctuation_ke_interval, save_interval
integer :: n_waves_in_r

! Functions of r:
real(8), allocatable, dimension(:) :: rho_midplane, Omega_K, ci, H

! To satisfy numerical vertical hydrostatic balance:
real(8), allocatable, dimension(:, :, :) :: pressure, dpdz
logical :: perturb, wavy_perturbation

logical :: apply_fargo_trick, apply_fargo_correction, integer_shifts = .false., &
     suppress_z_derivatives_when_nz_not_1_arg = .false., periodic_z_arg = .false., &
     use_supplied_dt = .false.

! To satisfy numerical centrifugal balance:
real(8), allocatable, dimension(:, :, :) :: p_r_space, dpdr_r_space, dpdr
real(8) :: uphi_squared

logical :: apply_artificial_pressure, get_lambda_max_ap
real(8) :: lambda_max_ap
logical :: use_rsize_for_domain, use_Manger_p
real(8) :: C_ap

! Viscosity stuff, e.g., LES, etc.
logical :: apply_viscosity
integer :: viscosity_type
real(8) :: C_Smag, nu_molecular, nu_b_molecular, Pr_molecular
real(8) :: C_DDSV ! Coefficient for dilatation-dependent shear (will be read from input file).

integer :: rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced

! Fluctuation kinetic energy:
real(8) :: afke_z_middle, afke_r_middle, afke_phi_middle, &
           afke_z_upper, afke_r_upper, afke_phi_upper, &
           afke_total_middle, afke_total_upper

logical :: name_using_step

logical, parameter :: stretched_r = .false., stretched_z = .false.
real(8), parameter :: r0_unif = 0.d0, z0 = 0.d0
integer, parameter :: nr_u = 0, nz_u = 0

character(3) :: sponge_type
logical :: apply_sponge
real(8) :: rho1, rho2, d1, d2, tau_decay
!integer :: n_decay_steps

character(25) :: filename_given
logical :: negative_density_flag

! Arguments for tecplot output calls.
integer :: iphi_plot, iz_plot
real(8) :: L_scale, T_scale, rho_scale

logical :: plot_pert = .false., plot_curl_rho_u = .true.

! For computing Reynolds phi averages from which Favre phi-time averaged means and stresses
! can be computed using TOOLS/phi_time_Favre_averages.f90.  This is valid only in a stationary
! state.
logical :: output_phi_Reynolds_averages
real(8) :: time_interval_for_phi_Reynolds_averages
real(8) :: next_t_for_phi_Reynolds_averages
logical, parameter :: average_in_z = .false. ! For Reynolds averages

namelist /vsi_3D_input/ &
     restart, file_version, isothermal, nr, nz, nphi, phi_max_over_pi, cfl, dt_min, nsteps, &
     tecplot_interval, output_profiles, profiles_interval, &
     output_phi_Reynolds_averages, time_interval_for_phi_Reynolds_averages, &
     fluctuation_ke_interval, save_interval, &
     perturb, wavy_perturbation, &
     zmax_over_H0, &
     use_rsize_for_domain, rsize_over_H0, &
     n_waves_in_r, use_Manger_p, &
     apply_pade_filter, eps_or_tau, eps_filter, &
     tau_filter, apply_artificial_pressure, C_ap, apply_viscosity, viscosity_type, &
     C_DDSV, name_using_step, &
     rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced, apply_fargo_trick, apply_fargo_correction, &
     apply_sponge, sponge_type, rho1, rho2, d1, d2, tau_decay

if (my_node .eq. 0) then
   print *, ' node 0: First executable of vsi_3D'
   print *, ' my_node = ', my_node, ' num_nodes = ', num_nodes   
end if

! Read namelist input:
if (my_node .eq. 0) then
   print *, ' node 0: about to open and read namelist file for app_vsi_3D'      
   open (unit = lun_general_purpose, file = 'input_file', form = 'formatted', &
        status = 'old')
   read (lun_general_purpose, nml = vsi_3D_input)
   close (lun_general_purpose)

   open(unit = lun_info, file = 'info.txt', form = 'formatted', &
        status = 'unknown', access = 'append')
   write(lun_info, *)' restart                   = ', restart
   write(lun_info, *)' file_version              = ', file_version
   write(lun_info, *)' isothermal                = ', isothermal
   write(lun_info, *)' nr, nz, nphi              = ', nr, nz, nphi
   write(lun_info, *)' phi_max_over_pi           = ', phi_max_over_pi
   write(lun_info, *)' cfl                       = ', cfl
   write(lun_info, *)' dt_min                    = ', dt_min
   write(lun_info, *)' nsteps                    = ', nsteps
   write(lun_info, *)' tecplot_interval          = ', tecplot_interval
   write(lun_info, *)' profiles_interval         = ', profiles_interval
   write(lun_info, *)' fluctuation_ke_interval   = ', fluctuation_ke_interval
   write(lun_info, *)' save_interval             = ', save_interval
   write(lun_info, *)''
   write(lun_info, *)' perturb                   = ', perturb
   write(lun_info, *)' wavy_perturbation         = ', wavy_perturbation
   write(lun_info, *)''
   write(lun_info, *)' zmax_over_H0              = ', zmax_over_H0
   write(lun_info, *)''
   write(lun_info, *)' use_rsize_for_domain      = ',  use_rsize_for_domain
   write(lun_info, *)' rsize_over_H0             = ', rsize_over_H0
   write(lun_info, *)' n_waves_in_r              = ', n_waves_in_r
   write(lun_info, *)' use_Manger_p              = ', use_Manger_p
   write(lun_info, *)''
   write(lun_info, *)' apply_pade_filter         = ', apply_pade_filter
   write(lun_info, *)' eps_or_tau                = ', eps_or_tau
   write(lun_info, *)' eps_filter                = ', eps_filter
   write(lun_info, *)' tau_filter                = ', tau_filter
   write(lun_info, *)''
   write(lun_info, *)' apply_artificial_pressure = ', apply_artificial_pressure
   write(lun_info, *)' C_ap                      = ', C_ap
   write(lun_info, *)''
   write(lun_info, *)' apply_viscosity           = ', apply_viscosity
   write(lun_info, *)' viscosity_type            = ', viscosity_type
   write(lun_info, *)' C_DDSV                    = ', C_DDSV
   write(lun_info, *)''
   write(lun_info, *)' name_using_step           = ', name_using_step
   write(lun_info, *)''
   write(lun_info, *)' rmin_BC                   = ', rmin_BC
   write(lun_info, *)' rmax_BC                   = ', rmax_BC
   write(lun_info, *)' zmin_BC                   = ', zmin_BC
   write(lun_info, *)' zmax_BC                   = ', zmax_BC
   write(lun_info, *)' ibalanced                 = ', ibalanced
   write(lun_info, *)''
   write(lun_info, *)' apply_fargo_trick         = ', apply_fargo_trick
   write(lun_info, *)' apply_fargo_correction    = ', apply_fargo_correction
   write(lun_info, *)''
   write(lun_info, *)' apply_sponge              = ', apply_sponge
   write(lun_info, *)' rho1                      = ', rho1
   write(lun_info, *)' rho2                      = ', rho2
   write(lun_info, *)' tau_decay                 = ', tau_decay
   close(lun_info)
end if

#ifdef mpi_code
   call mpi_bcast(restart,                   1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(file_version,              1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(isothermal,                1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(nr,                        1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nz,                        1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(nphi,                      1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(phi_max_over_pi,           1, mpi_double,  0, mpi_comm_world, ier)   
   call mpi_bcast(cfl,                       1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(dt_min,                    1, mpi_double,  0, mpi_comm_world, ier)   
   call mpi_bcast(nsteps,                    1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(tecplot_interval,          1, mpi_integer, 0, mpi_comm_world, ier)

   call mpi_bcast(output_profiles,           1, mpi_logical, 0, mpi_comm_world, ier)   
   call mpi_bcast(profiles_interval,         1, mpi_integer, 0, mpi_comm_world, ier)

   call mpi_bcast(output_phi_Reynolds_averages, 1, mpi_logical, 0, mpi_comm_world, ier)   
   call mpi_bcast(time_interval_for_phi_Reynolds_averages, 1, mpi_double, &
        0, mpi_comm_world, ier)   
   
   call mpi_bcast(fluctuation_ke_interval,              1, mpi_integer, 0, mpi_comm_world, ier)   
   call mpi_bcast(save_interval,             1, mpi_integer, 0, mpi_comm_world, ier)
   
   call mpi_bcast(perturb,                   1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(wavy_perturbation,         1, mpi_logical, 0, mpi_comm_world, ier)
   
   call mpi_bcast(zmax_over_H0,              1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(use_rsize_for_domain,      1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(rsize_over_H0,             1, mpi_double,  0, mpi_comm_world, ier)   
   call mpi_bcast(n_waves_in_r,              1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(use_Manger_p,              1, mpi_logical, 0, mpi_comm_world, ier)   

   call mpi_bcast(apply_pade_filter,         1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(eps_or_tau,                3, mpi_character, 0, mpi_comm_world, ier)
   call mpi_bcast(eps_filter,                1, mpi_double,  0, mpi_comm_world, ier)   
   call mpi_bcast(tau_filter,                1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(apply_artificial_pressure, 1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(C_ap,                      1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(apply_viscosity,           1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(viscosity_type,            1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(C_DDSV,                    1, mpi_double,  0, mpi_comm_world, ier)   
   call mpi_bcast(name_using_step,           1, mpi_logical, 0, mpi_comm_world, ier)

   call mpi_bcast(rmin_BC,                   1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(rmax_BC,                   1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(zmin_BC,                   1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(zmax_BC,                   1, mpi_integer, 0, mpi_comm_world, ier)
   call mpi_bcast(ibalanced,                 1, mpi_integer, 0, mpi_comm_world, ier)

   call mpi_bcast(apply_fargo_trick,         1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(apply_fargo_correction,    1, mpi_logical, 0, mpi_comm_world, ier)

   call mpi_bcast(apply_sponge,              1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(sponge_type,               3, mpi_character, 0, mpi_comm_world, ier)   
   call mpi_bcast(rho1,                      1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(rho2,                      1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(d1,                        1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(d2,                        1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(tau_decay,                 1, mpi_double,  0, mpi_comm_world, ier)
#endif


! The basic state follows Nelson et al. (2013).

! We can set three things to unity:
!GM    = 1.0d0
pi         = 4.0*ATAN(1.0d0)
T_orbital0 = 1.0d0
H0         = 1.0d0
rho0       = 1.0d0

! This is a parameter:
H0_over_r0 = 0.10d0
r0 = H0 / H0_over_r0

! Consequence of above:
Omega0 = 2.d0 * pi / T_orbital0
GM     = Omega0**2 * r0**3
V_Kep0 = SQRT(GM/r0)
c0     = H0_over_r0 * V_Kep0

! Gravity set-up parameters:
gravity_flag = 1
i_thin       = 0
i_no_radial  = 0

! Exponents:
if (use_Manger_p) then
   pexp = -2.d0/3.d0  ! for density.  Manger
else
   pexp = -3.d0/2.d0  ! for density Nelson
end if
qexp = -1.0d0  ! for temperature
c_exponent = 0.5d0*qexp
H_exponent = (qexp + 3.d0) / 2.d0

zmin = - zmax_over_H0 * H0
zmax =   zmax_over_H0 * H0

! Most amplified mode according to Matt:
lambda_r = pi * abs(qexp) * H0_over_r0 * H0

! Domain size in r:
if (.not. use_rsize_for_domain) then 
   rsize_over_H0 = n_waves_in_r * lambda_r
end if

rmin   = r0 - 0.5d0 * rsize_over_H0
rmax   = r0 + 0.5d0 * rsize_over_H0

phi_min = 0.0d0
phi_max = phi_max_over_pi * pi

! These calls should be in every application subroutine:
call set_up_domain_mesh_and_partition(rmin, rmax, zmin, zmax, phi_min, phi_max, nr, nz, nphi, &
     suppress_z_derivatives_when_nz_not_1_arg, periodic_z_arg, &
     stretched_r, stretched_z, r0_unif, nr_u, z0, nz_u)               
call set_up_thermal_parameters(gamma, isothermal)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: Returned from set_up_thermal_parameters'
   read(5, *)
#endif
call set_up_boundary_conditions(rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: Returned from set_up_boundary_conditions'
#endif   

call activate_gravity(gravity_flag, i_thin, i_no_radial, GM, .false., 0.0d0)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: Returned from activate_gravity'
#endif   
if (apply_fargo_trick) then
   ! The arguments are Boolean:
   call activate_fargo(integer_shifts, apply_fargo_correction)
   if (my_node .eq. 0) print *, ' node 0: Returned from activate_fargo'
end if

call activate_pade_filter(apply_pade_filter, eps_or_tau, eps_filter, &
     tau_filter, filter_relative_to_basic_state)
if (apply_artificial_pressure) then
   call activate_artificial_pressure (C_ap)
end if

if (apply_viscosity) then
   C_Smag     = 0.22d0**2
   nu_molecular   = 0.d0
   nu_b_molecular = 0.d0 ! bulk viscosity
   Pr_molecular = 0.d0
   ! subroutine activate_viscosity(apply_viscosity_arg, viscosity_type_arg, isothermal_arg, &
   ! nu_molecular_arg, nu_b_molecular_arg, Pr_molecular_arg, gamma_arg, C_DDSV_arg, C_Smag_arg)
   
   call activate_viscosity(apply_viscosity, viscosity_type, isothermal, &
     nu_molecular, nu_b_molecular, Pr_molecular, gamma, C_DDSV, C_Smag)
end if  

! Needed for non-reflective BC (for isothermal):
if (isothermal) then
   do iz = 1, nz
      d_ci_dr_inner(iz) = c0 * c_exponent / r0 * (rgrid(1 )/r0)**(c_exponent - 1.d0)
      d_ci_dr_outer(iz) = c0 * c_exponent / r0 * (rgrid(nr)/r0)**(c_exponent - 1.d0)
   end do
end if

! ------------------------------------------------------
! Compute basic state/initial condition and put it in q:
! ------------------------------------------------------
! For an non-restart run q will be the initial condition.
! For BOTH restart and non-restart runs, part of q will
! be used to store the basic state.

! We compute the basic state so it can be used to output fluctuations.
! Previously we used to store the basic state
allocate(rho_midplane(nr), Omega_K(nr), ci(nr), H(nr))
allocate(pressure(sr:er, sphi:ephi, nz), dpdz(sr:er, sphi:ephi, nz))

! For numerical centrifugal balance:
allocate(p_r_space   (sphi:ephi, sz_r:ez_r, nr))
allocate(dpdr_r_space(sphi:ephi, sz_r:ez_r, nr))
allocate(dpdr        (sr:er, sphi:ephi, nz))

! Functions of r for the basic state:
do ir = 1, nr
   rho_midplane(ir) = rho0 * (rgrid(ir) / r0)**pexp
   Omega_K     (ir) = SQRT(GM/rgrid(ir)**3)
   ci          (ir) = c0 * (rgrid(ir)/r0)**c_exponent
   H           (ir) = H0 * (rgrid(ir)/r0)**H_exponent
end do

! Needed for non-reflective BC:
do iz = 1, nz
   d_ci_dr_inner(iz) = c0 * c_exponent / r0 * (rgrid(1 )/r0)**(c_exponent - 1.d0)
   d_ci_dr_outer(iz) = c0 * c_exponent / r0 * (rgrid(nr)/r0)**(c_exponent - 1.d0)
end do

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         R_spherical = SQRT(rgrid(ir)**2 + zgrid(iz)**2)
         arg         = GM / ci(ir)**2 * (1.d0/R_spherical - 1.d0/rgrid(ir))
         rho_temporary          = rho_midplane(ir) * EXP(arg)
         pressure(ir, iphi, iz) = rho_temporary * ci(ir)**2
      end do
   end do
end do

! Initial condition:
if (.not. restart) then
   call pade_diff_z(mr*mphi, pressure, dpdz)

   ! Set rho to satisfy numerical vertical hydrostatic balance:
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            ! Required since gz = 0 at the midplane:
            if (gz(ir, iz) .ne. 0.d0) then
               q(ir, iphi, iz, irho) = dpdz(ir, iphi, iz)/gz(ir, iz)
            else
               q(ir, iphi, iz, irho) = rho_midplane(ir)
            end if
            ! Correct the pressure:
            pressure(ir, iphi, iz) = q(ir,iphi,iz,irho) * ci(ir)**2
         end do
      end do
   end do

   ! Set uphi to satisfy numerical centrifugal balance:
   call transpose_z_to_r (1, pressure, p_r_space)
   call pade_diff_bundle(mphi*mz_r, nr, Ji_r, p_r_space, dpdr_r_space)
   call transpose_r_to_z (1, dpdr_r_space, dpdr)   
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            q(ir, iphi, iz, rmom) = 0.d0
            q(ir, iphi, iz, zmom) = 0.0d0

            ! term1 = (p + qexp) * (H(ir)/rgrid(ir))**2
            ! term2 = 1.d0 + qexp
            ! R_spherical = SQRT(rgrid(ir)**2 + zgrid(iz)**2)
            ! term3 = - qexp*rgrid(ir) / R_spherical
            ! Omega = Omega_K(ir) * (term1 + term2 + term3)**(0.5d0)
            ! uphi  = Omega * rgrid(ir)

            ! Look at NRR-3 notes:
            uphi_squared = rgrid(ir)*(1.d0/q(ir,iphi,iz,irho)*dpdr(ir,iphi,iz) - gr(ir,iz))
            q(ir, iphi, iz, amom) = q(ir, iphi, iz, irho) * sqrt(uphi_squared) * rgrid(ir)

            ! Set energy in case we run adiabatically:
            ! q(ir, iphi, iz, ener) = pressure(ir,iphi,iz)/gm1 + 0.5d0 * q(ir, iphi, iz, irho) * uphi**2
            q(ir, iphi, iz, ener) = 0.d0
         end do
      end do
   end do
end if
! Finished with basic state

deallocate(pressure, dpdz, p_r_space, dpdr_r_space, dpdr)

! Store what is in q as the basic state:
call store_basic_state(filter_relative_to_basic_state)

if (restart) then
   call read_restart (file_version, istep_of_restart_file, t_of_restart_file)
else
   ! Assign initial condition using the basic state:
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            q(ir, iphi, iz, irho) = rho_basic(ir, iz)
            q(ir, iphi, iz, rmom) = 0.d0
            q(ir, iphi, iz, zmom) = 0.d0                        
            q(ir, iphi, iz, amom) = rho_basic(ir, iz) * uphi_basic(ir, iz) * rgrid(ir)

            ! Set energy in case we run adiabatically:
            ! q(ir, iphi, iz, ener) = pressure(ir,iphi,iz)/gm1
            q(ir, iphi, iz, ener) = 0.d0
         end do
      end do
   end do
end if
! Finished with initial condition.

if (apply_sponge) then
   call activate_sponge(sponge_type, rho1, rho2, rho_basic, d1, d2, tau_decay)   
end if

! Needed for the isothermal option.  ci_squared_initial sits in
! module thermal_parameters.
if (isothermal) then
   do iz = 1, nz
      do ir = 1, nr
         ci_squared_initial(ir, iz) = ci(ir)**2
      end do
   end do
end if

if (.not. restart) then
   istep0 = 0
   t      = 0.0d0
else
   istep0 = istep_of_restart_file
   t      = t_of_restart_file
end if

if (perturb) then
   if (wavy_perturbation) then
      if (my_node .eq. 0) print *, ' app_vsi_3D: Calling add_wavy_perturbation'
      call add_wavy_perturbation(c0, H0_over_r0, H0, lambda_r)
      if (my_node .eq. 0) print *, ' app_vsi_3D: Returned from wavy_perturbation'
   else
      call add_random_perturbation(c0)
      if (my_node .eq. 0) print *, ' app_vsi_3D: Returned from add_random_perturbation'      
   end if
end if

! These routines are currently only serial.
#ifndef mpi_code
   ! call output_gravity_profile_at_mid_radius
   ! call vertical_profiles(ir_mid, 1, 'midr', t) ! mid ir
   ! call radial_profiles  (iz_mid, 1, 'midz', t) ! midplane
   ! call radial_profiles  (1,      1, 'botz', t) ! bottom
   ! call radial_profiles  (nz,     1, 'topz', t) ! top
#endif

call sanity_check

if (my_node .eq. 0) print *, ' node 0: in vsi_3D Mb allocated = ', float(8*n_words_allocated)/1.d6

if (my_node .eq. 0) then
   print *, ' calling tecplot_fluc_vel_meridional_plane'
end if

iphi_plot    = 1
iz_plot      = iz_mid
!iz_plot      = 192
L_scale      = 1.0d0
T_scale      = 1.0d0
rho_scale    = 1.0d0
plot_curl_rho_u = .true.
call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, iz_plot, &
     L_scale, T_scale, rho_scale, t, istep0)
plot_curl_rho_u = .false.
call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, iz_plot, &
     L_scale, T_scale, rho_scale, t, istep0)
call tecplot_fluc_vel_meridional_plane         (q, iphi_plot,          L_scale, T_scale, t, istep0)
call tecplot_meridional_plane                  (q, iphi_plot,          L_scale, T_scale, t, istep0)
if (nphi .ne. 1) call tecplot_horizontal_plane (q,            iz_plot,          T_scale, t, istep0)
call phi_Favre_averaged_statistics(istep0, t)

if (output_phi_Reynolds_averages) then
   call phi_Reynolds_averages(istep0, t, average_in_z)
   next_t_for_phi_Reynolds_averages = t + time_interval_for_phi_Reynolds_averages
end if

call write_ave_fluctuation_ke(t)

! Time stepping loop:
do istep = istep0 + 1, istep0 + nsteps
   call rk4(rhs, q, cfl, t, dt, use_supplied_dt, negative_density_flag)
   call output_conservation_diagnostics(t, q)
   
   if (my_node .eq. 0) then
      print *, ' routine app_vsi_3D, node 0: finished istep = ', istep, ' t = ', t
   end if

   !call print_max(q)

   if (mod(istep - istep0, tecplot_interval) .eq. 0) then
      call tecplot_meridional_plane (q, 1, 1.0d0, 1.0d0, t, istep)
      if (nphi .ne. 1) call tecplot_horizontal_plane (q, iz_plot, 1.0d0, t, istep)

      plot_curl_rho_u = .true.      
      call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, iz_plot, &
           L_scale, T_scale, rho_scale, t, istep)
      plot_curl_rho_u = .false.            
      call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, iz_plot, &
           L_scale, T_scale, rho_scale, t, istep)      
      call tecplot_fluc_vel_meridional_plane(q, 1, 1.0d0, 1.0d0, t, istep)
      
      ! Debug:
      ! Plot the artificial pressure and dilatation in the plane in which the eigenvalue lambda_max_ap
      ! is the largest:
      !if (apply_artificial_pressure) then
      !   get_lambda_max_ap = .true.
      !   if (my_node .eq. 0) print *, ' calling get_artificial_pressure' 
      !   call get_artificial_pressure(get_lambda_max_ap, q, lambda_max_ap)
      !   if (my_node .eq. 0) print *, ' calling tecplot_p_art_in_meridional_plane'
      !   call tecplot_scalar_in_meridional_plane('p_art', p_art, iphi_max_ap, 1.0d0, 1.0d0, &
      !        1.0d0, t, istep)
      !   call tecplot_scalar_in_meridional_plane('dilat', dil,   iphi_max_ap, 1.0d0, 1.0d0, &
      !                                            1.0d0, t, istep)         
      !end if
   end if

   if ((output_phi_Reynolds_averages) .and. (t .ge. next_t_for_phi_Reynolds_averages)) then
      call phi_Reynolds_averages(istep, t, average_in_z)
      next_t_for_phi_Reynolds_averages = t + time_interval_for_phi_Reynolds_averages
   end if   

   if (mod(istep - istep0, fluctuation_ke_interval) .eq. 0) then
      call write_ave_fluctuation_ke(t)
   end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node = 0, vsi_3D: about to check for a save'
#endif
   
   if (mod(istep - istep0, save_interval) .eq. 0) then
      call write_save_file(istep, t, filename_given)
   end if

   ! Detect too small time step:
   if (istep .gt. istep0+1) then
      if (dt .lt. dt_min) then
         if (my_node .eq. 0) then
            print *, ' my_node = ', my_node
            print *, ' dt < dt_min'
            print *, ' dt          = ', dt
            print *, ' dt_min      = ', dt_min
            print *, ' dt_previous = ', dt_previous
         end if
      end if
   end if

   if ((dt .lt. dt_min) .or. (negative_density_flag)) then
      call tecplot_meridional_plane (q, 1, 1.0d0, 1.0d0, t, istep)
      if (nphi .ne. 1) call tecplot_horizontal_plane (q, iz_plot, 1.0d0, t, istep)

      call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, iz_plot, &
           L_scale, T_scale, rho_scale, t, istep)      
      call tecplot_fluc_vel_meridional_plane(q, 1, 1.0d0, 1.0d0, t, istep)
      close(lun_history(1))
      ! 1 = istatus : Abnormal return.
      call terminate_with_save(1, istep, t)
   end if
   
   dt_previous = dt
end do

if (my_node .eq. 0) then
   print *, ' ***** average stepping time per step = ', total_cpu_for_stepping/nsteps
end if

close(lun_history(1))
! 0 = istatus: Normal return
! Need the minus 1 since a fortran do loop increments the counter at the end
call terminate_with_save(0, istep-1, t)

end subroutine vsi_3D

!----------------------------------------------------------------------------------85

subroutine add_wavy_perturbation(c0, H0_over_r0, H0, lambda_r)

use q_array
use grid
use dof_indices
use partition_data
#ifdef mpi_code
   use mpi, only: mpi_comm_world
#endif

! ifort compatibility library so we can use the same random number generator as gfortran.   
#ifdef ifort
   use IFPORT 
#endif
   
implicit none
real(8), intent(in) :: c0, H0_over_r0, H0, lambda_r
integer, parameter :: nkr = 11, nkz = 12, nkphi_dim = 12
integer :: nkphi
real(8) :: kr(nkr), kz(nkz), kphi(nkphi_dim)
! Local:
integer :: ir, iz, iphi, ikr, ikz, ikphi, ier
real(8) :: pi, eps, lambda_z

real(8) :: lambda_r_modulation, kr_modulation, lambda_z_modulation, kz_modulation, &
     r_modulation, z_modulation, modulation
real(8) :: rho, uz, ur, uphi
real(8) :: uz_prime, ur_prime, uphi_prime
real(8) :: phase_ur, phase_uz, phase_uphi
complex(8) :: exp_arg, exp_fac, ic = CMPLX(0.0d0, 1.0d0)

complex(8), allocatable, dimension(:,:,:) :: A_ur, A_uz, A_uphi

if (my_node .eq. 0) print *, ' modified add_perturbation_3D has been called'

if (nphi .eq. 1) then
   nkphi = 1
else
   nkphi = nkphi_dim
end if

allocate(A_ur(nkphi,nkz,nkr), A_uz(nkphi,nkz,nkr), A_uphi(nkphi,nkz,nkr))

! Establish seed for random number generator different for each processor:
call srand(987654 + my_node*769)
pi = 4.0d0 * atan(1.0d0)
eps = 1.d-3*c0

! The purpose of the modulation is to make normal velocities = 0 at boundaries.
! The modulation is a half sin.
lambda_r_modulation = 2.d0 * (rgrid(nr) - rgrid(1))
kr_modulation       = 2.d0 * pi / lambda_r_modulation
lambda_z_modulation = 2.d0 * (zgrid(nz) - zgrid(1))
kz_modulation       = 2.d0 * pi / lambda_z_modulation

! Note: nkr = 11
kr(1) = 2.d0 * pi / lambda_r
do ikr = 2, 8
   kr(ikr) = ikr * kr(1)
end do
! Sub-harmonics:
kr(9)  = kr(1)/2.d0
kr(10) = kr(1)/3.d0
kr(11) = kr(1)/4.d0

lambda_z = zgrid(nz) - zgrid(1)
kz(1) = 2.d0 * pi / lambda_z
do ikz = 2, nkz
   kz(ikz) = ikz * kz(1)
end do

do ikphi = 1, nkphi
   kphi(ikphi) = ikphi - 1
end do

if (my_node .eq. 0) then
   print *, ' lambda_r = ', lambda_r
   print *, ' # of grid points per lambda_r = ', nr * lambda_r / (rmax - rmin)
   print *, ' # of waves in r domain = ', (rgrid(nr) - rgrid(1)) / lambda_r
end if

#ifdef mpi_code
   call mpi_barrier(mpi_comm_world, ier)
#endif

! This entire section in incorrect for what you want to do.


! Get complex amplitudes for each mode:
do ikr = 1, nkr
   do ikz = 1, nkz
      do ikphi = 1, nkphi
         ! Zero argument means next random number in the sequence.
         phase_ur   = rand(0) * 2.d0 * pi                  
         phase_uz   = rand(0) * 2.d0 * pi
         phase_uphi = rand(0) * 2.d0 * pi

         A_uz  (ikphi, ikz, ikr) = cos(phase_uz)   + ic*sin(phase_uz)
         A_ur  (ikphi, ikz, ikr) = cos(phase_ur)   + ic*sin(phase_ur)
         A_uphi(ikphi, ikz, ikr) = cos(phase_uphi) + ic*sin(phase_uphi)
      end do
   end do
end do
   
do iz = 1, nz
   ! print *, ' in point loop iz = ', iz
   ! print *, ' sphi = ', sphi, ' ephi = ', ephi
   z_modulation = sin(kz_modulation*(zgrid(iz) - zgrid(1)))
   do iphi = sphi, ephi
      do ir = sr, er
         ! r_modulation = sin(kr_modulation*(rgrid(ir) - rgrid(1)))
         ! modulation = r_modulation * z_modulation
         modulation = z_modulation
         
         rho  = q(ir, iphi, iz, irho)
         uz   = q(ir, iphi, iz, zmom) / rho
         ur   = q(ir, iphi, iz, rmom) / rho
         uphi = q(ir, iphi, iz, amom) / rho / rgrid(ir)

         ! Sum over modes:
         uz_prime   = 0.0d0
         ur_prime   = 0.0d0
         uphi_prime = 0.0d0
         do ikr = 1, nkr
            do ikz = 1, nkz
               do ikphi = 1, nkphi
                  exp_arg = ic * (kr(ikr)*rgrid(ir) + kz(ikz)*zgrid(iz) + &
                       kphi(ikphi)*phi_grid(iphi))
                  exp_fac = EXP(exp_arg)
                  uz_prime   = uz_prime   + A_uz  (ikphi,ikz,ikr) * exp_fac
                  ur_prime   = ur_prime   + A_ur  (ikphi,ikz,ikr) * exp_fac
                  uphi_prime = uphi_prime + A_uphi(ikphi,ikz,ikr) * exp_fac
                  ! print *, ' phase_uz = ', phase_uz
                  ! print *, ' A_uz = ', A_uz, ' pi = ', pi
               end do
            end do
         end do ! sum over modes
         uz_prime   = eps * modulation * uz_prime
         ur_prime   = eps * modulation * ur_prime
         uphi_prime = eps * modulation * uphi_prime

         q(ir, iphi, iz, zmom) = rho*(uz   + uz_prime)
         q(ir, iphi, iz, rmom) = rho*(ur   + ur_prime)
         q(ir, iphi, iz, amom) = rho*(uphi + uphi_prime)*rgrid(ir)
      end do
   end do
end do

end subroutine add_wavy_perturbation

!----------------------------------------------------------------------------------85

subroutine add_random_perturbation(c0)

use q_array
use grid
use dof_indices
use partition_data
use math_constants
#ifdef mpi_code
   use mpi, only: mpi_comm_world
#endif

! ifort compatibility library so we can use the same random number generator as gfortran.   
#ifdef ifort
   use IFPORT 
#endif
   
implicit none
real(8), intent(in) :: c0
integer :: ir, iz, iphi, ier
real(8) :: eps

real(8) :: lambda_r_modulation, kr_modulation, lambda_z_modulation, kz_modulation, &
     r_modulation, z_modulation, modulation
real(8) :: rho, uz, ur, uphi
real(8) :: uz_prime, ur_prime, uphi_prime

if (my_node .eq. 0) print *, ' random_perturbation has been called'

! Establish seed for random number generator:
call srand(987654 + 789*my_node)
eps = 1.d-3*c0

! The purpose of the modulation is to make normal velocities = 0 at boundaries.
! The modulation is a half sin.
lambda_r_modulation = 2.d0 * (rgrid(nr) - rgrid(1))
kr_modulation       = 2.d0 * pi / lambda_r_modulation
lambda_z_modulation = 2.d0 * (zgrid(nz) - zgrid(1))
kz_modulation       = 2.d0 * pi / lambda_z_modulation

#ifdef mpi_code
   call mpi_barrier(mpi_comm_world, ier)
#endif

do iz = 1, nz, 2
   ! print *, ' in point loop iz = ', iz
   ! print *, ' sphi = ', sphi, ' ephi = ', ephi
   z_modulation = sin(kz_modulation*(zgrid(iz) - zgrid(1)))
   do iphi = sphi, ephi, 2
      do ir = sr, er, 2
         r_modulation = sin(kr_modulation*(rgrid(ir) - rgrid(1)))
         modulation = r_modulation * z_modulation
         modulation = z_modulation
         modulation = 1.d0
         
         rho  = q(ir, iphi, iz, irho)
         uz   = q(ir, iphi, iz, zmom) / rho
         ur   = q(ir, iphi, iz, rmom) / rho
         uphi = q(ir, iphi, iz, amom) / rho / rgrid(ir)

         uz_prime   = eps*rand(0)*modulation
         ur_prime   = eps*rand(0)*modulation
         uphi_prime = eps*rand(0)*modulation

         q(ir, iphi, iz, zmom) = rho*(uz   + uz_prime)
         q(ir, iphi, iz, rmom) = rho*(ur   + ur_prime)
         q(ir, iphi, iz, amom) = rho*(uphi + uphi_prime)*rgrid(ir)
      end do
   end do
end do

end subroutine add_random_perturbation

!----------------------------------------------------------------------------------85





