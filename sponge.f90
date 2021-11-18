!----------------------------------------------------------------------------------85

! Form of sponge: - (q - q_basic) / tau_sponge

!----------------------------------------------------------------------------------85

module sponge
   logical :: apply_sponge
   ! This determines tau_sponge:   tau_sponge = n_decay_steps * dt
   !integer :: n_decay_steps

   real(8) :: tau_decay
   ! This is 1 where the sponge is to be applied
   !         0 where not applied
   !         and varies smoothly in between using Melander's function.
   real(8), allocatable, dimension(:, :) :: sponge_mask
end module sponge

!----------------------------------------------------------------------------------85

! I put on rho_basic in the argument list so that the user knows that this routine
! should be called after the basic state has been set.

subroutine activate_sponge(sponge_type, rho1, rho2, rho_basic, d1, d2, tau_decay_arg)

! Sponge where the basic state density is smaller than rho1.
! rho1 and rho2 define a smooth transition.
!
! For sponge_type = 'bsd'   (basic state density)
! Sponge amplitude =
! 1 for rho_basic < rho1
! transition for rho1 < rho_basic < rho2
! 0 for rho_basic > rho2

! For sponge_type = 'drw'   (distance from radial wall)
! Sponge amplitude =
! 1 for radial distance from radial wall < d1
! transition for d1 < radial distance from radial wall < d2
! 0 for radial distance from radial wall > d2

! For sponge_type = 'daw'   (distance from all walls)
! Sponge amplitude =
! 1 for normal distance from wall < d1
! transition for d1 < radial distance from wall < d2
! 0 for normal distance from wall > d2

! 
! For sponge_type = 'cor' (sponge at left bottom and left top corner)
! Sponge amplitude = 1
! 1 for normal distance from wall < d1
! transition for d1 < radial distance from wall < d2
! 0 for normal distance from wall > d2


use sponge
use grid
use partition_data
implicit none
! The sponge is fully applied for rho < rho1
! Not applied at all for rho > rho2
! And, transitions in between.
character(3),                  intent(in) :: sponge_type
real(8),                       intent(in) :: rho1, rho2
real(8), dimension(sr:er, nz), intent(in) :: rho_basic
real(8),                       intent(in) :: d1, d2
real(8),                       intent(in) :: tau_decay_arg

! Local:
integer :: iz, ir
real(8) :: Melander3 ! called function
real(8) :: distance_from_left, distance_from_right, &
     distance_from_bottom, distance_from_top, distance

! The 1 is to prevent conflict with variables in module grid:
real(8) :: dr1, dz1

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' First executable in activate sponge'
   print *, ' tau_decay_arg = ', tau_decay_arg
end if
#endif

apply_sponge  = .true.
!n_decay_steps = n_decay_steps_arg
tau_decay = tau_decay_arg

allocate(sponge_mask(sr:er, nz))

#ifdef debug_print
if (my_node .eq. 0) then
   print *, ' Finished allocating sponge mask'
end if
#endif

if (sponge_type .eq. 'bsd') then
   do iz = 1, nz
      do ir = sr, er
         sponge_mask(ir,iz) = Melander3(rho_basic(ir,iz), rho1, rho2)
         if (sponge_mask(ir,iz) .lt. 0.d0) then
            print *, ' sponge_mask < 0'
         end if
      end do
   end do
else if (sponge_type .eq. 'drw') then
   do iz = 1, nz
      do ir = sr, er
         distance_from_left  = abs(rgrid(ir) - rgrid(1 ))
         distance_from_right = abs(rgrid(ir) - rgrid(nr))
         distance = min(distance_from_left, distance_from_right)
         sponge_mask(ir,iz) = Melander3(distance, d1, d2)
         if (sponge_mask(ir,iz) .lt. 0.d0) then
            print *, ' sponge_mask < 0'
         end if
      end do
   end do
else if (sponge_type .eq. 'daw') then
   do iz = 1, nz
      do ir = sr, er
         distance_from_left   = abs(rgrid(ir) - rgrid(1 ))
         distance_from_right  = abs(rgrid(ir) - rgrid(nr))
         distance_from_bottom = abs(zgrid(iz) - zgrid(1 ))
         distance_from_top    = abs(zgrid(iz) - zgrid(nz))
         
         distance = min(distance_from_left,   distance_from_right, &
                        distance_from_bottom, distance_from_top)
         
         sponge_mask(ir,iz) = Melander3(distance, d1, d2)
         if (sponge_mask(ir,iz) .lt. 0.d0) then
            print *, ' sponge_mask < 0'
         end if
      end do
   end do
else if (sponge_type .eq. 'cor') then
   ! Top left and bottom left sponge:
   do iz = 1, nz
      do ir = sr, er
         dr1 = rgrid(ir) - rgrid(1)
         ! Bottom:
         dz1 = zgrid(iz) - zgrid(1)
         distance_from_bottom = sqrt(dr1**2 + dz1**2)

         ! Top
         dz = zgrid(iz) - zgrid(nz)         
         distance_from_top = sqrt(dr1**2 + dz1**2)
         
         distance = min(distance_from_bottom, distance_from_top)
         
         sponge_mask(ir,iz) = Melander3(distance, d1, d2)
         if (sponge_mask(ir,iz) .lt. 0.d0) then
            print *, ' sponge_mask < 0'
         end if
      end do
   end do   
else
   print *, 'Unknown sponge type'
   call terminate_with_no_save(1)
end if

call tecplot_axisymmetric_scalar('SpongeMask', sponge_mask, &
     1.0d0, 1.0d0, 1.0d0, 0.d0)

end subroutine activate_sponge

!----------------------------------------------------------------------------------85

subroutine add_sponge_term(dt, q, qdot)

use grid
use dof_indices, only: irho, zmom, amom, rmom, ener
use partition_data
use sponge
use basic_state
use thermal_parameters
implicit none
real(8), intent(in) :: dt
real(8), intent(in),    dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8), intent(inout), dimension(sr:er, sphi:ephi, nz, ndof) :: qdot

! Local:
integer :: iz, iphi, ir
real(8) :: tau_inv

#ifdef debug_print
   if (my_node .eq. 0) print *, ' First executable in add_sponge_term'
#endif

tau_inv = 1.d0 / tau_decay
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         qdot(ir,iphi,iz,irho) = qdot(ir,iphi,iz,irho)-&
              sponge_mask(ir,iz)*(q(ir,iphi,iz,irho) - rho_basic(ir,iz))*tau_inv
         qdot(ir,iphi,iz,rmom) = qdot(ir,iphi,iz,rmom)-&
              sponge_mask(ir,iz)*(q(ir,iphi,iz,rmom)                 )*tau_inv
         qdot(ir,iphi,iz,zmom) = qdot(ir,iphi,iz,zmom)-&
              sponge_mask(ir,iz)*(q(ir,iphi,iz,zmom)                 )*tau_inv
         qdot(ir,iphi,iz,amom) = qdot(ir,iphi,iz,amom)-&
              sponge_mask(ir,iz)*(q(ir,iphi,iz,amom)-amom_basic(ir,iz))*tau_inv
      end do
   end do
end do

if (.not. isothermal) then
   do iz = 1, nz
      do iphi = sphi, ephi
         do ir = sr, er
            !qdot(ir,iphi,iz,ener) = qdot(ir,iphi,iz,ener)-sponge_mask(ir,iz)*(q(ir,iphi,iz,ener)                 )               *tau_inv
         end do
      end do
   end do
end if

end subroutine add_sponge_term

!----------------------------------------------------------------------------------85

