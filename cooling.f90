!----------------------------------------------------------------------------------85

subroutine cooling_terms (q, p, qdot)

use physical_constants
use dof_indices
use thermal_parameters
use grid
use partition_data
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(in)    :: q    ! field
real(8), dimension(sr:er, sphi:ephi, nz),       intent(in)    :: p    ! pressure
real(8), dimension(sr:er, sphi:ephi, nz, ndof), intent(inout) :: qdot 

! Locals:
real(8), dimension(sr:er, sphi:ephi, nz) :: T
real(8), dimension(sr:er, sphi:ephi, nz) :: capNH

capNH = 0.0d0

do iphi = sphi, ephi
   do ir = sr, er
      capNH(ir, iphi, 1) 
         rho = q(ir, iphi, iz, irho) ! density
         T(ir, iphi, iz) = p(ir, iphi, iz) / (rho * Rgas)
         nH = rho / mH
         capNH = ??

         
          Lambda_CO  = 9.79d-12 * nH/capNH * T**5.5d0
            Lambda_H2O = 4.84d-11 * nH/capNH * T**6.0d0
         end do
      end do
   end do
end subroutine cooling_terms

!----------------------------------------------------------------------------------85

!----------------------------------------------------------------------------------85

