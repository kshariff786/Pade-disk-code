!----------------------------------------------------------------------------------85

module activate_routines

contains

! Activation routines for:
! gravity, artificial_pressure, pade_filter, fargo, plotting_shift

!----------------------------------------------------------------------------------85

subroutine activate_fargo(integer_shifts_arg, apply_fargo_extra_operator_arg)

! A call to this routine should be done in the set-up portion of the user's app.
! It activates the Fargo method for reducing the time step if it is dominated by
! Keplerian advection due to a central mass.

! We use the suffix _arg so that the subroutine arguments don't clash with
! variables in module fargo_trick_or_plotting_shift.

! integer_shifts_arg: If .true., we will use the original scheme of Masset which
! employs integer shifts which can undergo a jump w.r.t. r.  If .false., we will
! apply real valued shifts.

! apply_fargo_extra_operator_arg:  If .true. we will apply the extra operator that
! arises from the chain-rule.  If .false., we will not apply this operator.

! GM : Product of Newton's gravitational constant and the mass of the central
! object.

use grid
use partition_data
use fargo_or_plotting_shift
use control_parameters
use transposes_of_q_and_qdot
use math_constants
use gravity, only: gravity_flag, GM
implicit none
logical :: integer_shifts_arg, apply_fargo_extra_operator_arg

! We should not apply Fargo if nphi = 1:
if (nphi .eq. 1) then
   if (my_node .eq. 0) then
      print *, ' Message from subroutine activate_fargo:'
      print *, ' Cannot apply the Fargo trick when nphi = 1'
      print *, ' Returning without activating Fargo'
   end if
   return
end if

apply_fargo                = .true.
integer_shifts             = integer_shifts_arg   
apply_fargo_extra_operator = apply_fargo_extra_operator_arg   

allocate(Omega_fargo(nr), uphi_fargo_subtract(nr), fargo_factor(nr))
allocate(nshift(nr), inew_phi(nr, nphi))
n_words_allocated = n_words_allocated + 3*nr + nr/2 + nr*nphi/2

! The second test is there because set_up_plotting_shift could have allocated it.
if ( (.not. integer_shifts) .and. (.not. allocated(complex_shift_factor)) ) then
   allocate(complex_shift_factor(nphi/2+1, sr:er))
   n_words_allocated = n_words_allocated + 2*(nphi/2+1)*mr
end if

if (.not. have_Omega_bar) call setup_Omega_bar_for_fargo(gravity_flag, GM)

if (.not. fft_has_been_initialized) then
#ifdef fftw
   call make_plans_for_fftw (q_phi_space)
   fft_has_been_initialized = .true.
#else
   call initialize_rogallo_fft_for_fargo
   fft_has_been_initialized = .true.
#endif
end if

twopi_over_Delta_phi_domain = 2.d0 * pi / Delta_phi_domain

end subroutine activate_fargo

!----------------------------------------------------------------------------------85

subroutine activate_plotting_shift(Omega0_plotting_shift_arg, t_of_previous_shift_arg)

! A call to this routine is optional.
! Set-up for Fargo trick for phi advection.

use grid
use partition_data
use fargo_or_plotting_shift
use control_parameters
use transposes_of_q_and_qdot
use math_constants
implicit none
logical :: apply_plotting_shift_arg
real(8) :: Omega0_plotting_shift_arg, t_of_previous_shift_arg

apply_plotting_shift  = .true.
Omega0_plotting_shift = Omega0_plotting_shift_arg
t_of_previous_shift   = t_of_previous_shift_arg

if (.not. allocated(complex_shift_factor) ) then
   allocate(complex_shift_factor(nphi/2+1, sr:er))
   n_words_allocated = n_words_allocated + 2*(nphi/2+1)*mr
end if

if (.not. fft_has_been_initialized) then
#ifdef fftw
   call make_plans_for_fftw (q_phi_space)
   fft_has_been_initialized = .true.
#else
   call initialize_rogallo_fft_for_fargo
   fft_has_been_initialized = .true.
#endif
end if

twopi_over_Delta_phi_domain = 2.d0 * pi / Delta_phi_domain

!print *, ' In set_up_plotting_shift'
!print *, ' twopi_over_Delta_phi_domain = ', twopi_over_Delta_phi_domain
!print *, ' enter anything to continue'
!read (5, *)

end subroutine activate_plotting_shift

!----------------------------------------------------------------------------------85

subroutine activate_gravity(gravity_flag_in, ithin, i_no_radial, GM_in, &
     finite_sphere_arg, R_sphere_arg)

! This is an optional routine.  The defult is no gravity.
! This routine is called to activate or deactivate gravity via
! gravity_flag.

! ithin: Flag for thin disk gravity.
! i_no_radial : Flag to suppress radial gravity for the thin disk option.
!               This is used to test vertical hydrostatic balance without
!               the radial momentum equation.
! GM_in : Value of GM.

! Thin disk gravity is: g = -GM/r^2 * (rhat + (z/r) zhat) where rhat and zhat are
! unit vectors.

! g = - grad Phi_g where Phi_g is the gravitational potential energy.

use grid
use gravity
use control_parameters, only: n_words_allocated
use partition_data
use math_constants
use fargo_or_plotting_shift
implicit none
integer :: gravity_flag_in, ithin, i_no_radial
logical :: finite_sphere_arg
real(8) :: GM_in, R_sphere_arg

! Local:
integer iz, ir
real(8) :: R_spherical, g_mag, cos_theta, sin_theta
real(8) :: rho_sphere, vol, GM_of_R

#ifdef debug_print
   if (my_node .eq. 0) print *, ' gravity_set_up has been called with GM_in = ', GM_in
#endif

finite_sphere = finite_sphere_arg
R_sphere      = R_sphere_arg   

gravity_flag = gravity_flag_in

allocate (gr(nr, nz), gz(nr, nz), Phi_g(nr, nz))
n_words_allocated = n_words_allocated + 3*nr*nz

GM = GM_in

! We do this because we would like to suggest Fargo to the user
! if we determine that it might be advantageous.
if(.not. have_Omega_bar) then
   ! Note: If igravity = 0 then this routine sets Omega_bar to zero so it can
   ! be used by the time-step determination routine.
   call setup_Omega_bar_for_fargo(gravity_flag, GM)
end if


if (gravity_flag .ne. 1) then
   ! Zero gravity:
   do iz = 1, nz
      do ir = 1, nr
         Phi_g(ir, iz) = 0.0d0
         gr   (ir, iz) = 0.0d0
         gz   (ir, iz) = 0.0d0
      end do
   end do
   return
end if

if (finite_sphere) then
   rho_sphere = GM / (4.d0/3.d0 * pi * R_sphere**3)
   do iz = 1, nz
      do ir = 1, nr
         R_spherical = SQRT(rgrid(ir)**2 + zgrid(iz)**2)
         if (R_spherical .lt. R_sphere) then
            vol = 4.d0/3.d0 * pi * R_spherical
            GM_of_R = rho_sphere * vol
         else
            GM_of_R = GM
         end if
         g_mag       = GM_of_R/ R_spherical**2
         cos_theta   = rgrid(ir) / R_spherical
         sin_theta   = zgrid(iz) / R_spherical
         ! We are choosing the sign so this equals the gravitational potential energy
         ! such that g = - grad Phi_g.
         Phi_g(ir, iz) = - GM_of_R/R_spherical
         gr   (ir, iz) = - g_mag * cos_theta
         gz   (ir, iz) = - g_mag * sin_theta
      end do
   end do
   return
end if
         
if ((ithin .eq. 1) .and. (i_no_radial .eq. 1)) then
   do iz = 1, nz
      do ir = 1, nr
         gr(ir, iz) = 0.0d0
         gz(ir, iz) = -GM/rgrid(ir)**2 * (zgrid(iz) / rgrid(ir))

         ! This is for a thin disk:
         R_spherical = rgrid(ir) * (1.d0 + 0.5d0*zgrid(iz)**2/rgrid(ir)**2)
         Phi_g (ir, iz) = - GM/R_spherical
      end do
   end do
else if ((ithin .eq. 1) .and. (i_no_radial .ne. 1)) then
   do iz = 1, nz
      do ir = 1, nr
         gr(ir, iz) = -GM/rgrid(ir)**2
         gz(ir, iz) = gr(ir, iz) * (zgrid(iz) / rgrid(ir))
         
         ! This is for a thin disk:
         R_spherical = rgrid(ir) * (1.d0 + 0.5d0*zgrid(iz)**2/rgrid(ir)**2)
         Phi_g (ir, iz) = - GM/R_spherical
      end do
   end do
else ! ithin = 0
   do iz = 1, nz
      do ir = 1, nr
         R_spherical = SQRT(rgrid(ir)**2 + zgrid(iz)**2)
         g_mag       = GM/ R_spherical**2
         cos_theta   = rgrid(ir) / R_spherical
         sin_theta   = zgrid(iz) / R_spherical
         ! We are choosing the sign so this equals the gravitational potential energy
         ! such that g = - grad Phi_g.
         Phi_g(ir, iz) = - GM/R_spherical
         gr   (ir, iz) = - g_mag * cos_theta
         gz   (ir, iz) = - g_mag * sin_theta
      end do
   end do
end if

! Debug:
!open (unit = 1, file = 'gz.dat', form = 'formatted', status = 'unknown')
!do iz = 1, nz
!   write (1, "(2(1x, e12.5))") zgrid(iz), gz(1, iz)
!end do
!close (1)

allocate (u_Kepler(nr), Omega_Kepler(nr))
n_words_allocated = n_words_allocated + 2*nr

if (gravity_flag .eq. 1) then
   do ir = 1, nr
      Omega_Kepler(ir) = SQRT((GM / rgrid(ir)**3))
      u_Kepler    (ir) = Omega_Kepler(ir) * rgrid(ir)
   end do
end if

end subroutine activate_gravity

!----------------------------------------------------------------------------------85

subroutine activate_artificial_pressure (C_ap_arg)

use partition_data
use artificial_pressure_module
use grid
use control_parameters

real(8), intent(in) :: C_ap_arg  ! constant for artificial pressure term

! Locals:
integer :: ir, iphi, iz 

apply_artificial_pressure = .true.
C_ap = C_ap_arg

if (.not. allocated(dil)) then
   allocate(dil(sr:er, sphi:ephi, nz))
   n_words_allocated = n_words_allocated + mr*mphi*nz
end if

! To do: Do you need all these?  Yes you do, in order to take derivatives along
! the various directions.
allocate(p_art          (sr:er, sphi:ephi, nz))
allocate(p_art_r_space  (sphi:ephi, sz_r:ez_r, nr))
allocate(p_art_phi_space(sr:er, sz_phi:ez_phi, nphi))
n_words_allocated = n_words_allocated + 3*mr*mz_phi*nphi

if (.not. partitioned) then
   print *, ' You need to partition the grid before calling activate_artificial_pressure'
   call terminate_with_no_save(1)
end if

if (.not. gridded) then
   print *, ' You need to generate the grid before calling activate_artificial_pressure'
   call terminate_with_no_save(1)
end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: About to return from activate artificial pressure'
#endif

end subroutine activate_artificial_pressure

!----------------------------------------------------------------------------------85

subroutine activate_pade_filter(apply_pade_filter_arg, eps_or_tau_arg, &
     eps_filter_arg, tau_filter_arg, filter_relative_to_basic_state_arg)

use partition_data
use pade_filter
use logical_units
implicit none
logical,      intent(in) :: apply_pade_filter_arg
character(3), intent(in) :: eps_or_tau_arg
real(8),      intent(in) :: eps_filter_arg, tau_filter_arg
logical,      intent(in) :: filter_relative_to_basic_state_arg

apply_pade_filter = apply_pade_filter_arg
eps_or_tau = eps_or_tau_arg
eps_filter = eps_filter_arg
tau_filter = tau_filter_arg
filter_relative_to_basic_state = filter_relative_to_basic_state_arg

if (my_node .eq. 0) then
   open(unit = lun_info, file = 'info.txt', form = 'formatted', status = 'unknown', &
        access = 'append')
   write (lun_info, 1) eps_or_tau, eps_filter, tau_filter
1  format (' In activate_pade_filter eps_or_tau = ', a3, ' eps_filter = ', e12.5, ' tau_filter = ', e12.5)
   close(lun_info)
end if

end subroutine activate_pade_filter

!----------------------------------------------------------------------------------85

end module activate_routines

!----------------------------------------------------------------------------------85
