!----------------------------------------------------------------------------------85

subroutine add_gravity_terms(q, qdot)

! This assumes axisymmetric (r, z) gravity and adds the r and z components of
! gravity terms to the momentum and energy equation.

! If we are running with isothermal = true, the energy equation is not needed
! but its presence here does no harm.

use grid, only: ndof, nz, nphi, nr
use dof_indices, only: irho, zmom, amom, rmom, ener
use gravity, only: gr, gz
use partition_data
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q, qdot

! Local:
integer :: iz, iphi, ir
real(8) :: fr, fz, ur, uz

do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         ! Gravitational forces per unit volume:
         fr = q(ir, iphi, iz, irho) * gr(ir, iz)
         fz = q(ir, iphi, iz, irho) * gz(ir, iz)
         ur = q(ir, iphi, iz, rmom) / q(ir, iphi, iz, irho)
         uz = q(ir, iphi, iz, zmom) / q(ir, iphi, iz, irho)
         qdot(ir,iphi,iz,rmom) = qdot(ir,iphi,iz,rmom) + fr 
         qdot(ir,iphi,iz,zmom) = qdot(ir,iphi,iz,zmom) + fz
         qdot(ir,iphi,iz,ener) = qdot(ir,iphi,iz,ener) + fr*ur + fz*uz
      end do
   end do
end do

end subroutine add_gravity_terms

!----------------------------------------------------------------------------------85

subroutine output_gravity_profile_at_mid_radius

! This is a serial routine.

use grid
use gravity
use logical_units
implicit none

! Local:
integer :: iz

open (unit = lun_profile(1), file = 'gz_gr_vs_z.dat', form = 'formatted', &
     status = 'unknown')
do iz = 1, nz
   write (lun_profile(1), 1) zgrid(iz), gz(ir_mid, iz), gr(ir_mid, iz)
   1 format (3(1x, e12.5))
end do

close (lun_profile(1))
end

!----------------------------------------------------------------------------------85
