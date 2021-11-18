!----------------------------------------------------------------------------------85

subroutine app_zvi

use dof_indices, only: irho, amom, zmom, rmom, ener
use grid, only: nr, nz, nphi, rmin, rmax, zmin, zmax, phi_min, phi_max, &
     ir_mid, iz_mid, rgrid, zgrid, Ji_r
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
logical :: apply_pade_filter, filter_relative_to_basic_state = .true.
character(3) :: eps_or_tau
real(8) :: eps_filter, tau_filter

! For thermal set-up:
real(8) :: gamma = 5.d0/3.d0
logical :: isothermal

! For gravity set-up
real(8) :: GM
integer :: gravity_flag, i_thin, i_no_radial

real(8) :: pexp, qexp, c_exponent, H_exponent ! exponents for rho, temp, etc.
! Formula temps.:
real(8) :: Omega, rho_temporary, term1, term2, term3, R_spherical, arg

integer :: iz, ir, iphi, ier
logical :: output_profiles, output_profiles_now
integer :: tecplot_interval, profiles_interval, rms_interval, save_interval
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
real(8) :: C_ap

integer :: rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced
real(8) :: uz_rms, ur_rms, uphi_rms

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

real(8), dimension(5) :: var_abs_max

namelist /zvi_input/ &
     restart, file_version, nr, nz, nphi, phi_max_over_pi, cfl, dt_min, nsteps, &
     tecplot_interval, output_profiles, profiles_interval, &
     output_phi_Reynolds_averages, time_interval_for_phi_Reynolds_averages, &
     rms_interval, save_interval, perturb, &
     zmax_over_H0, rsize_over_H0, &
     apply_pade_filter, eps_or_tau, eps_filter, &
     tau_filter, apply_artificial_pressure, C_ap, &
     rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced, apply_fargo_trick, apply_fargo_correction, &
     apply_sponge, sponge_type, rho1, rho2, d1, d2, tau_decay

if (my_node .eq. 0) then
   print *, ' node 0: First executable of zvi'
   print *, ' my_node = ', my_node, ' num_nodes = ', num_nodes   
end if

#ifdef mpi_code
   call mpi_barrier(mpi_comm_world, ier)
#endif

! Read namelist input:
if (my_node .eq. 0) then
   print *, ' node 0: about to open and read namelist file for app_zvi'      
   open (unit = lun_general_purpose, file = 'input_file', form = 'formatted', status = 'old')
   read (lun_general_purpose, nml = zvi_input)
   close (lun_general_purpose)

   print *, ' restart                   = ', restart
   print *, ' file_version              = ', file_version
   print *, ' nr, nz, nphi              = ', nr, nz, nphi
   print *, ' phi_max_over_pi           = ', phi_max_over_pi
   print *, ' cfl                       = ', cfl
   print *, ' dt_min                    = ', dt_min
   print *, ' nsteps                    = ', nsteps
   print *, ' tecplot_interval          = ', tecplot_interval
   print *, ' profiles_interval         = ', profiles_interval
   print *, ' rms_interval              = ', rms_interval
   print *, ' save_interval             = ', save_interval
   print *, ''
   print *, ' perturb                   = ', perturb
   print *, ' zmax_over_H0              = ', zmax_over_H0
   print *, ' rsize_over_H0             = ', rsize_over_H0
   print *, ''
   print *, ' apply_pade_filter         = ', apply_pade_filter
   print *, ' eps_or_tau                = ', eps_or_tau
   print *, ' eps_filter                = ', eps_filter
   print *, ' tau_filter                = ', tau_filter
   print *, ''
   print *, ' apply_artificial_pressure = ', apply_artificial_pressure
   print *, ' C_ap                      = ', C_ap
   print *, ''

   print *, ''
   print *, ' rmin_BC                   = ', rmin_BC
   print *, ' rmax_BC                   = ', rmax_BC
   print *, ' zmin_BC                   = ', zmin_BC
   print *, ' zmax_BC                   = ', zmax_BC
   print *, ' ibalanced                 = ', ibalanced
   print *, ''
   print *, ' apply_fargo_trick         = ', apply_fargo_trick
   print *, ' apply_fargo_correction    = ', apply_fargo_correction
   print *, ''
   print *, ' apply_sponge              = ', apply_sponge
   print *, ' rho1                      = ', rho1
   print *, ' rho2                      = ', rho2
end if

#ifdef mpi_code
   call mpi_bcast(restart,                   1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(file_version,              1, mpi_integer, 0, mpi_comm_world, ier)
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
   
   call mpi_bcast(rms_interval,              1, mpi_integer, 0, mpi_comm_world, ier)   
   call mpi_bcast(save_interval,             1, mpi_integer, 0, mpi_comm_world, ier)
   
   call mpi_bcast(perturb,                   1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(zmax_over_H0,              1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(rsize_over_H0,             1, mpi_double,  0, mpi_comm_world, ier)   

   call mpi_bcast(apply_pade_filter,         1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(eps_or_tau,                3, mpi_character, 0, mpi_comm_world, ier)
   call mpi_bcast(eps_filter,                1, mpi_double,  0, mpi_comm_world, ier)   
   call mpi_bcast(tau_filter,                1, mpi_double,  0, mpi_comm_world, ier)
   call mpi_bcast(apply_artificial_pressure, 1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(C_ap,                      1, mpi_double,  0, mpi_comm_world, ier)

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

if (my_node .eq. 0) then
   print *, ' routine app_zvi, my_node = ', my_node, ' tau_filter = ', tau_filter
end if

! We are running adiabatically:
isothermal = .false.

! The basic state follows Nelson et al. (2013).

! We can set three things to unity:
!GM    = 1.0d0
pi         = 4.0*ATAN(1.0d0)
T_orbital0 = 1.0d0
H0         = 1.0d0
rho0       = 1.0d0

! This is a parameter:
H0_over_r0 = 0.05d0
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
! For ZVI there is no radial temperature gradient.
pexp = 0.d0  ! for density
qexp = 0.d0  ! for temperature
c_exponent = 0.5d0*qexp
H_exponent = (qexp + 3.d0) / 2.d0

zmin = - zmax_over_H0 * H0
zmax =   zmax_over_H0 * H0

rmin   = r0 - 0.5d0 * rsize_over_H0
rmax   = r0 + 0.5d0 * rsize_over_H0

phi_min = 0.0d0
phi_max = phi_max_over_pi * pi

! These calls should be in every application subroutine:
call set_up_domain_mesh_and_partition(rmin, rmax, zmin, zmax, phi_min, phi_max, nr, nz, nphi, &
     suppress_z_derivatives_when_nz_not_1_arg, periodic_z_arg, &
     stretched_r, stretched_z, r0_unif, nr_u, z0, nz_u)               
call set_up_thermal_parameters(gamma, isothermal)
call set_up_boundary_conditions(rmin_BC, rmax_BC, zmin_BC, zmax_BC, ibalanced)

call activate_gravity(gravity_flag, i_thin, i_no_radial, GM, .false., 0.0d0)
if (apply_fargo_trick) then
   ! The first two arguments are Boolean:
   call activate_fargo(integer_shifts, apply_fargo_correction)
   if (my_node .eq. 0) print *, ' node 0: Returned from activate_fargo'
end if

call activate_pade_filter(apply_pade_filter, eps_or_tau, eps_filter, &
     tau_filter, filter_relative_to_basic_state)
if (apply_artificial_pressure) then
   call activate_artificial_pressure (C_ap)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' returned from activate_artificial_pressure'
#endif
end if

if (restart) then
   call read_restart (file_version, istep_of_restart_file, t_of_restart_file)
end if

allocate(rho_midplane(nr), Omega_K(nr), ci(nr), H(nr))
allocate(pressure(sr:er, sphi:ephi, nz), dpdz(sr:er, sphi:ephi, nz))

! For numerical centrifugal balance:
allocate(p_r_space   (sphi:ephi, sz_r:ez_r, nr))
allocate(dpdr_r_space(sphi:ephi, sz_r:ez_r, nr))
allocate(dpdr        (sr:er, sphi:ephi, nz))

! Functions of r for the basic state:
do ir = 1, nr
   rho_midplane(ir) = rho0
   Omega_K     (ir) = SQRT(GM/rgrid(ir)**3)
   ci          (ir) = c0
   H           (ir) = H0 * (rgrid(ir)/r0)**H_exponent
end do

#ifdef debug_print
   if (my_node .eq. 0) print *, ' Finished computing functions of r for basic state'
#endif

! Needed for non-reflective BC:
!do iz = 1, nz
!   d_ci_dr_inner(iz) = c0 * c_exponent / r0 * (rgrid(1 )/r0)**(c_exponent - 1.d0)
!   d_ci_dr_outer(iz) = c0 * c_exponent / r0 * (rgrid(nr)/r0)**(c_exponent - 1.d0)
!end do

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
#ifdef debug_print
   if (my_node .eq. 0) print *, ' In app_zvi computing dp_dz for initial condition'
#endif
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

            q(ir, iphi, iz, ener) = 0.d0
         end do
      end do
   end do

   if (.not. isothermal) then
      do iz = 1, nz
         do iphi = sphi, ephi
            do ir = sr, er
               q(ir, iphi, iz, ener) = pressure(ir,iphi,iz)/gm1            
            end do
         end do
      end do
   end if
end if
! Finished with initial condition.

deallocate(pressure, dpdz, p_r_space, dpdr_r_space, dpdr)

if (.not. restart) then
   ! Store the initial condition as the basic state:
   call store_basic_state(filter_relative_to_basic_state)
end if

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
   call zvi_add_random_perturbation(c0)
   if (my_node .eq. 0) print *, ' app_zvi: Returned from zvi_add_random_perturbation'      
end if

call sanity_check

if (my_node .eq. 0) print *, ' node 0: in app_zvi Mb allocated = ', float(8*n_words_allocated)/1.d6

! --------------
! Initial plots:
! --------------
iphi_plot    = 1
iz_plot      = iz_mid
L_scale      = 1.0d0
T_scale      = 1.0d0
rho_scale    = 1.0d0

call tecplot_meridional_plane(q, iphi_plot, L_scale, T_scale, t, istep0)

if (my_node .eq. 0) then
   print *, ' calling tecplot_fluc_vel_meridional_plane'
end if

call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, &
     iphi_plot, iz_plot, L_scale, T_scale, rho_scale, t, istep0)

#ifdef debug_print
   if (my_node .eq. 0) print *, ' rank = 0: app_zvi: tecplot_fluc_vel_meridional plane'
#endif

call tecplot_fluc_vel_meridional_plane         (q, iphi_plot,          L_scale, T_scale, t, istep0)

#ifdef debug_print
   if (my_node .eq. 0) print *, ' rank = 0: Returned from tecplot_fluc_vel_meridional plane'
#endif

#ifdef debug_print
   if (my_node .eq. 0) print *, ' rank = 0: Returned from tecplot_meridional plane'
#endif

if (nphi .ne. 1) then
   call tecplot_horizontal_plane (q,            iz_plot,          T_scale, t, istep0)
end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' rank = 0: About to call output_rms'
#endif

call output_rms(t, uz_rms, ur_rms, uphi_rms, output_profiles)

#ifdef debug_print
   if (my_node .eq. 0) print *, ' rank = 0: About to call phi_Favre_averaged_statistics'
#endif
   
call phi_Favre_averaged_statistics(istep0, t)

if (output_phi_Reynolds_averages) then
   call phi_Reynolds_averages(istep0, t, average_in_z)
   next_t_for_phi_Reynolds_averages = t + time_interval_for_phi_Reynolds_averages
end if

if (my_node .eq. 0) then
   open(unit = lun_history(1), file = 'uz_ur_uphi_rms.his', form = 'formatted', &
        status = 'unknown', access = 'append')   
   write(lun_history(1), "(4(1x, e12.5))") t, uz_rms, ur_rms, uphi_rms
   ! In case the code bombs:
   close(lun_history(1))
   open(unit = lun_history(2), file = 'dt.his', form = 'formatted', &
        status = 'unknown', access = 'append')   
   close(lun_history(2))   
end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' rank = 0: About to begin time-stepping loop'
#endif

! Time stepping loop:
do istep = istep0 + 1, istep0 + nsteps
   call rk4(rhs, q, cfl, t, dt, use_supplied_dt, negative_density_flag)

   ! call output_conservation_diagnostics(t, q)
   
   if (my_node .eq. 0) then
      print *, ' routine app_vsi_3D, node 0: finished istep = ', istep, ' t = ', t
   end if

   call print_max(q, t, var_abs_max)

   if (my_node .eq. 0) then
      open(unit = lun_history(2), file= 'ur_max.dat', form = 'formatted', status = 'unknown', &
           access = 'append')
      write(lun_history(2), "(2(1x, e12.5))") t, var_abs_max(2)
      close(lun_history(2))
   end if

   if (mod(istep - istep0, tecplot_interval) .eq. 0) then
      call tecplot_meridional_plane (q, 1, 1.0d0, 1.0d0, t, istep)
      if (nphi .ne. 1) call tecplot_horizontal_plane (q, iz_plot, 1.0d0, t, istep)

      call tecplot_vort_and_dil_in_merid_horiz_planes(plot_pert, plot_curl_rho_u, q, iphi_plot, iz_plot, &
           L_scale, T_scale, rho_scale, t, istep)
      call tecplot_fluc_vel_meridional_plane(q, 1, 1.0d0, 1.0d0, t, istep)
   end if

   if ((output_phi_Reynolds_averages) .and. (t .ge. next_t_for_phi_Reynolds_averages)) then
      call phi_Reynolds_averages(istep, t, average_in_z)
      next_t_for_phi_Reynolds_averages = t + time_interval_for_phi_Reynolds_averages
   end if   

   if (output_profiles .and. (mod(istep - istep0, profiles_interval) .eq. 0)) then
      call output_rms(t, uz_rms, ur_rms, uphi_rms, output_profiles)
   end if

   if (mod(istep - istep0, rms_interval) .eq. 0) then
      output_profiles_now = .false.
      call output_rms(t, uz_rms, ur_rms, uphi_rms, output_profiles_now)      
      if (my_node .eq. 0) then         
         open(unit = lun_history(1), file = 'uz_ur_uphi_rms.his', form = 'formatted', &
              status = 'unknown', access = 'append')   
         write(lun_history(1), "(4(1x, e12.5))") t, uz_rms, ur_rms, uphi_rms
         ! In case the code bombs:
         close(lun_history(1))
         open(unit = lun_history(2), file = 'dt.his', form = 'formatted', &
              status = 'unknown', access = 'append')   
         write(lun_history(2), "(2(1x, e12.5))") t, dt
         ! In case the code bombs:
         close(lun_history(2))                  
      end if
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

! Need the minus 1 since a fortran do loop increments the counter at the end
close(lun_history(1))
! 0 = istatus: Normal return
call terminate_with_save(0, istep-1, t)

end subroutine app_zvi

!----------------------------------------------------------------------------------85

subroutine zvi_add_random_perturbation(c0)

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
call srand(987654 + 123*my_node)
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
   do iphi = sphi, ephi
      do ir = sr, er
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

end subroutine zvi_add_random_perturbation

!----------------------------------------------------------------------------------85





