!----------------------------------------------------------------------------------85

subroutine attack_dilatation(q, nsteps, cfl)

use grid
use partition_data
use artificial_pressure_module
use logical_units
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
integer, intent(in) :: nsteps
real(8), intent(in) :: cfl

! Local:
external rhs_attack_dilatation
integer :: istep
logical :: use_supplied_dt, negative_density
real(8) :: t, dt, dil_max_abs
integer :: ir_max, iz_max, iphi_max

if (my_node .eq. 0) then
   print *, ' In subroutine attack_dilatation'
   open(unit = lun_dil_abs_max, file = 'dil_abs_max.dat', form = 'formatted', &
        status = 'unknown')
end if

use_supplied_dt = .false.
t = 0.d0
istep = 0
call dilatation(q)
call max_abs_of_array(dil, dil_max_abs, ir_max, iz_max, iphi_max)
if (my_node .eq. 0) then
   print *, ' istep = ', istep, ' |dil|_max = ', dil_max_abs, &
        ' at (ir, iz, iphi) = ', ir_max, iz_max, iphi_max
   write(lun_dil_abs_max, "(i7, e12.5)") istep, dil_max_abs   
end if

do istep = 1, nsteps
   call rk4(rhs_attack_dilatation, q, cfl, t, dt, use_supplied_dt, negative_density)
   call dilatation(q)
   call max_abs_of_array(dil, dil_max_abs, ir_max, iz_max, iphi_max)
   if (my_node .eq. 0) then
      print *, ' istep = ', istep, ' |dil|_max = ', dil_max_abs, &
           ' at (ir, iz, iphi) = ', ir_max, iz_max, iphi_max
      write(lun_dil_abs_max, "(i7, e12.5)") istep, dil_max_abs         
   end if
end do
if (my_node .eq. 0) close(lun_dil_abs_max)

end subroutine attack_dilatation

!----------------------------------------------------------------------------------85

subroutine rhs_attack_dilatation(q, qdot, first_substep, t, t0, lambda_max_ap, &
     lambda_max_extra_operator)

! RHS of ODE system with only application of artificial pressure to drive down
! dilatation. t and t0 are needed to compute fargo_factor.

use grid
use partition_data
use dof_indices
use logical_units, only: lun_lambda
use math_constants
#ifdef mpi_code
   use mpi
#endif
use transposes_of_q_and_qdot
use fargo_or_plotting_shift
use artificial_pressure_module

implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q, qdot
logical, intent(in)                            :: first_substep
real(8), intent(in)                            :: t          ! needed to compute fargo_factor
real(8), intent(in )                           :: t0         ! needed to compute fargo_factor
real(8), intent(out)                           :: lambda_max_ap ! maximum eigenvalue/pi from all the terms
                                                             ! in the first sub-step.
real(8), intent(out)                           :: lambda_max_extra_operator

! Local:
! Make these global later:
real(8), dimension(sr:er,     sphi:ephi,     nz  ) :: dP_z_space
real(8), dimension(sphi:ephi, sz_r:ez_r,     nr  ) :: dP_r_space
real(8), dimension(sr:er,     sz_phi:ez_phi, nphi) :: dP_phi_space

integer :: iz, iphi, ir, ier

if (.not. apply_artificial_pressure) return
have_dilatation = .false.   
qdot = 0.0d0

if (first_substep) then
   add_fargo_extra_operator_now = .false.
else
   ! Note: apply_fargo_this_step will have been determined in the first sub-step.
   if (apply_fargo_this_step) then
      add_fargo_extra_operator_now = apply_fargo_extra_operator
      do ir = 1, nr
         fargo_factor(ir) = - (t - t0) * d_Omega_bar_dr(ir) / rgrid(ir)
      end do
   end if
end if

call get_artificial_pressure(first_substep, q, lambda_max_ap)

if ((my_node .eq. 0) .and. first_substep) then
   ! The items in the last two lines are in module artificial_pressure
   print *, ' In rhs_art_press_only'
   print *, ' lambda_max_ap = ', lambda_max_ap, &
            ' at ir = ', ir_max_ap, ' iz = ', iz_max_ap, ' iphi = ', iphi_max_ap, &
            ' beta_Delta_max = ', beta_Delta_max
   ! read(5, *)
end if

!!$if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
!!$   call pade_diff_z(nbundle_z, p_art, dP_z_space)
!!$   qdot(:, :, :, zmom) = qdot(:, :, :, zmom) - dP_z_space(:, :, :)   
!!$end if
!!$
!!$if (nr .ne. 1) then
!!$   call transpose_z_to_r(1, p_art, p_art_r_space)
!!$   call transpose_z_to_r(ndof, qdot, qdot_r_space)   
!!$   call pade_diff_bundle(mphi*mz_r, nr, Ji_r, p_art_r_space, dP_r_space)
!!$
!!$   do ir = 1, nr
!!$      do iz = sz_r, ez_r
!!$         do iphi = sphi, ephi
!!$            qdot_r_space(iphi,iz,rmom,ir) = qdot_r_space(iphi,iz,rmom,ir) - dP_r_space(iphi,iz,ir)
!!$         end do
!!$      end do
!!$   end do
!!$   ! Transpose qdot back to z space:
!!$   call transpose_r_to_z(ndof, qdot_r_space, qdot)   
!!$end if

if (nphi .ne. 1) then
   call transpose_z_to_phi(1, p_art, p_art_phi_space)
   call transpose_z_to_phi(ndof, qdot, qdot_phi_space)      
   call pade_diff_periodic(mr*mz_phi, nphi, dphi, p_art_phi_space, dP_phi_space)

   do iphi = 1, nphi
      do iz = sz_phi, ez_phi
         do ir = sr, er
            qdot_phi_space(ir,iz,amom,iphi) = qdot_phi_space(ir,iz,amom,iphi)-dP_phi_space(ir,iz,iphi)
         end do
      end do
   end do

   if (apply_fargo_this_step) then
      do iphi = 1, nphi
         do iz = sz_phi, ez_phi
            do ir = sr, er
               qdot_phi_space(ir,iz,rmom,iphi) = qdot_phi_space(ir,iz,rmom,iphi) - &
                    dP_phi_space(ir,iz,iphi) * fargo_factor(ir)
            end do
         end do
      end do
   end if
   
   call transpose_phi_to_z(ndof, qdot_phi_space, qdot)
end if

! The dilatation, etc. will now pertain to the middle of a substep.
have_dilatation               = .false.
lambda_max_extra_operator     = 0.d0 ! It is going to be small anyway.

! Sync. for safety:
#ifdef mpi_code
call mpi_barrier(mpi_comm_world, ier)
#endif

end subroutine rhs_attack_dilatation

!----------------------------------------------------------------------------------85
