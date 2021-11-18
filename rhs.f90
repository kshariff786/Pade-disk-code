!----------------------------------------------------------------------------------85

subroutine rhs(q, qdot, first_substep, t, t0, lambda_max, lambda_max_extra_operator)

! RHS of ODE system for governing equations for mass, momentum, and internal energy.
! t and t0 are needed to compute fargo_factor (if needed) which has (t - t0).

use grid
use partition_data
use dof_indices
use gravity, only: gravity_flag
use boundary_conditions
use logical_units, only: lun_lambda
use math_constants
use sponge
use rotating_frame
#ifdef mpi_code
   use mpi
#endif
use viscous ! Laminar terms and LES model
use transposes_of_q_and_qdot
use fargo_or_plotting_shift
! Location of maximum eigenvalue lambda_max_ap_global and bulk viscosity at
! this location.
use artificial_pressure_module, only: apply_artificial_pressure, have_dilatation, &
     ir_max_ap, iz_max_ap, iphi_max_ap, beta_Delta_max
use pade_filter

implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q, qdot
logical, intent(in)                            :: first_substep
real(8), intent(in)                            :: t          ! needed to compute fargo_factor which has (t-t0)
real(8), intent(in )                           :: t0         ! needed to compute fargo_factor
real(8), intent(out)                           :: lambda_max ! maximum eigenvalue/pi from all the terms
                                                             ! in the first sub-step.
real(8), intent(out)                           :: lambda_max_extra_operator

! Local:

! Eigenvalues for time step determination.
real(8) :: lambda_max_viscous, lambda_max_non_euler
real(8) :: lambda_max_ap_global, &
           lambda_max_extra_operator_global, lambda_max_viscous_global
real(8) :: lambda_max_euler, lambda_max_euler_fargo, lambda_max_euler_non_fargo
real(8) :: lambda_max_without_fargo, lambda_max_with_fargo

integer :: iz, iphi, ir, ier
logical :: worth_doing_fargo, dominated_by_phi

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' rhs has been called'
      print *, ' nz = ', nz
      print *, ' suppress_z_derivatives_when_nz_not_1 = ', suppress_z_derivatives_when_nz_not_1
   end if
#endif

lambda_max_extra_operator        = 0.d0
lambda_max_extra_operator_global = 0.d0
lambda_max_euler_fargo           = 0.d0
lambda_max_euler_non_fargo       = 0.d0

! These are used to possibly avoid calculating the dilatation in an expensive way.
! We are beginning a new substep so all old quantities should be recomputed if needed.   
have_dilatation               = .false.   
have_velocity_gradient_tensor = .false.
have_strain_tensor            = .false.

qdot = 0.0d0

if (first_substep) then
   if (apply_fargo) then
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In rhs, about to call get_lambda_euler_fargo'
#endif
call get_lambda_euler_fargo(q, lambda_max_euler_non_fargo, lambda_max_euler_fargo)
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In rhs, returned from get_lambda_euler_fargo'
#endif
   else
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In rhs, about to call get_lambda_euler_non_fargo'
#endif   
   call get_lambda_euler_non_fargo(q, lambda_max_euler_non_fargo, dominated_by_phi)
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In rhs, returned from get_lambda_euler_non_fargo'
#endif
      lambda_max_euler_fargo = 0.d0
   end if
   add_fargo_extra_operator_now = .false.
else
   ! Note: apply_fargo_this_step will have been determined in the first sub-step.
   if (apply_fargo_this_step) then
      add_fargo_extra_operator_now = apply_fargo_extra_operator
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In rhs, doing loop for fargo_factor'
#endif      
      do ir = 1, nr
         fargo_factor(ir) = - (t - t0) * d_Omega_bar_dr(ir) / rgrid(ir)
      end do
   end if
#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: In rhs, finished loop for fargo_factor'
#endif   
end if

! We do the viscous terms first because if we compute the strain rate or velocity gradient
! tensor it can be used to compute the dilatation for the artificial pressure.
if (apply_viscosity) then
   ! This is located in viscous_terms.f90
   call add_viscous_terms_to_qdot(t, q, qdot, lambda_max_viscous) ! In viscous_terms.f90
#ifdef debug_print
   if (my_node .eq. 0) print *, ' routine rhs, node 0: Returned from add_viscous_terms_to_qdot'
#endif
   if (first_substep) then
#ifdef mpi_code
      ! Note: Only proc. 0 will have lambda_max_viscous_global
      call mpi_reduce(lambda_max_viscous, lambda_max_viscous_global, 1, mpi_double, mpi_max, &
           0, mpi_comm_world, ier)
#else
      lambda_max_viscous_global = lambda_max_viscous
#endif
   end if
else
   lambda_max_viscous_global = 0.d0
end if

if (apply_viscosity .and. (my_node .eq. 0) .and. first_substep) then
   print *, ' lambda_max_viscous_global = ', lambda_max_viscous_global
end if

! --------------------------Artificial pressure ----------------------
if (apply_artificial_pressure) then
   ! This subroutine is in artificial_pressure.f90
   ! The actual adding of the artificial pressure will happen in the add_derivative routines.
#ifdef debug_print
   if (my_node .eq. 0) print *, ' routine rhs, node 0: Calling get_artificial_pressure'
#endif
   ! Adding in of artificial pressure happens later. This routine is in artificial_pressure.f90,
   ! It is only in the first-substep that we obtain lambda_max_ap_global.
   call get_artificial_pressure (first_substep, q, lambda_max_ap_global)
   if ((my_node .eq. 0) .and. first_substep) then
      ! The items in the last two lines are in module artificial_pressure
      print *, ' lambda_max_ap_global = ', lambda_max_ap_global, &
           ' at ir = ', ir_max_ap, ' iz = ', iz_max_ap, ' iphi = ', iphi_max_ap, &
           ' beta_Delta_max = ', beta_Delta_max
   end if
else
   lambda_max_ap_global = 0.d0
end if

if (first_substep) then
   if (my_node .eq. 0) then
      lambda_max_non_euler = max(lambda_max_viscous_global, lambda_max_ap_global)
      lambda_max_without_fargo = max(lambda_max_euler_non_fargo, lambda_max_non_euler)
      lambda_max_with_fargo    = max(lambda_max_euler_fargo,     lambda_max_non_euler)

      worth_doing_fargo = (lambda_max_with_fargo .lt. lambda_max_without_fargo) .and. &
           dominated_by_phi

      ! See if we should suggest Fargo to the user:
      if ((nphi .ne. 1) .and. (.not. apply_fargo)) then
         if (worth_doing_fargo) then
            print *, ' subroutine rhs: We suggest that you use Fargo'
            print *, ' lambda_max_with_fargo    = ', lambda_max_with_fargo
            print *, ' lambda_max_without_fargo = ', lambda_max_without_fargo

            print *, ' lambda_max_euler_non_fargo = ', lambda_max_euler_non_fargo
            print *, ' lambda_max_euler_fargo     = ', lambda_max_euler_fargo            
            print *, ' lambda_max_non_euler       = ', lambda_max_non_euler
         end if
      end if

      ! See if we should do Fargo this step:
      if (apply_fargo .and. worth_doing_fargo) then
         apply_fargo_this_step = .true.
         lambda_max = lambda_max_with_fargo
      else
         apply_fargo_this_step = .false.
         lambda_max = lambda_max_without_fargo
      end if
   end if ! my_node = 0
#ifdef mpi_code
   call mpi_bcast(apply_fargo_this_step, 1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(lambda_max,            1, mpi_double,  0, mpi_comm_world, ier)
#endif
end if ! first_substep

if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
#ifdef debug_print
   if (my_node .eq. 0) print *, ' routine rhs, node 0: Calling add_z_derivatives'
#endif   
   call add_z_derivatives(q, qdot)
end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: about to do if check on gravity flag in subroutine rhs'
#endif

! Gravity terms in the momentum equations and work term in the energy equation.   
if (gravity_flag .eq. 1) call add_gravity_terms(q, qdot) ! In gravity.f90

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: in rhs about to add centrifugal term'
#endif

! Centrifugal term in radial momentum eq.:
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         qdot(ir,iphi,iz,rmom) = qdot(ir,iphi,iz,rmom) + &
                  q(ir,iphi,iz,amom)**2/q(ir,iphi,iz,irho)/rgrid(ir)**3
      end do
   end do
end do

! Rotating frame terms if required:
if (apply_rotating_frame) call add_coriolis_and_centrifugal_terms(q, qdot)

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: in rhs finished adding centrifugal term'
#endif

if (nr .ne. 1) then
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: in rhs about to call add_radial_derivatives'
#endif      
   call add_radial_derivatives(q, qdot)
end if

if (nphi .ne. 1) then
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: about to call add_phi_derivatives'
#endif
   call add_phi_derivatives(q, qdot, lambda_max_extra_operator)
#ifdef mpi_code
   ! This reduction takes places to proc. 0:
   call mpi_reduce(lambda_max_extra_operator, lambda_max_extra_operator_global, 1, mpi_double, mpi_max, &
        0, mpi_comm_world, ier)
#else
   lambda_max_extra_operator_global = lambda_max_extra_operator
#endif
end if

! The dilatation, etc. will now pertain to the middle of a substep.
have_dilatation               = .false.   
have_velocity_gradient_tensor = .false.
have_strain_tensor            = .false.

! Sync. for safety:
#ifdef mpi_code
call mpi_barrier(mpi_comm_world, ier)
#endif

end subroutine rhs

!----------------------------------------------------------------------------------85
