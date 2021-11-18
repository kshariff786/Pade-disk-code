!----------------------------------------------------------------------------------85

subroutine rhs_shearing_box(q, qdot, first_substep, t, t0, lambda_max)

! RHS of ODE system for governing equations for mass, momentum, and internal energy
! for the shearing box.  Currently no viscous terms.

use grid
use partition_data
use dof_indices
use gravity, only: gravity_flag
use boundary_conditions
use logical_units, only: lun_lambda
use math_constants
#ifdef mpi_code
   use mpi
#endif   
use transposes_of_q_and_qdot
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

! Local:

integer :: iz, iphi, ir, ier

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' rhs has been called'
      print *, ' nz = ', nz
      print *, ' suppress_z_derivatives_when_nz_not_1 = ', suppress_z_derivatives_when_nz_not_1
   end if
#endif

qdot = 0.0d0

if (first_substep) then
   call get_lambda_euler_non_fargo(q, lambda_max_euler, dominated_by_phi)
#ifdef mpi_code
   call mpi_bcast(apply_fargo_this_step, 1, mpi_logical, 0, mpi_comm_world, ier)
   call mpi_bcast(lambda_max,            1, mpi_double,  0, mpi_comm_world, ier)
#endif   
end if

if (first_substep) then

      end if
   end if ! my_node = 0

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
have_grad_rho                 = .false.

! Sync. for safety:
#ifdef mpi_code
call mpi_barrier(mpi_comm_world, ier)
#endif

end subroutine rhs

!----------------------------------------------------------------------------------85
