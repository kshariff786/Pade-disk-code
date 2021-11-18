!----------------------------------------------------------------------------------85

module mass_diffusivity
   logical :: apply_mass_diffusivity, have_grad_rho
   integer :: mass_diffusivity_type
   real(8) :: C_DDMD ! coefficient for dilatation dependent mass diffusivity.
   real(8), allocatable, dimension(:, :, :) :: grad_rho_z, grad_rho_r, grad_rho_phi
   real(8), allocatable, dimension(:,:,:,:) :: mass_fluxes
   real(8), allocatable, dimension(:,:,:  ) :: mass_diff_term
end module mass_diffusivity

module mass_diffusivity_types
   integer, parameter :: Yoshizawa_mass_diffusivity = 1, dilatation_dependent_mass_diffusivity = 2
end module mass_diffusivity_types

!----------------------------------------------------------------------------------85

subroutine activate_mass_diffusivity(mass_diffusivity_type_arg, C_DDMD_arg)

use grid
use partition_data
use mass_diffusivity
use mass_diffusivity_types
use artificial_pressure_module, only: dil
use logical_units
use control_parameters, only: n_words_allocated
implicit none
integer :: mass_diffusivity_type_arg
real(8) :: C_DDMD_arg

if (my_node .eq. 0) then
   print *, ' activate_mass_diffusivity has been called:'
   write(lun_info, *) ' activate_mass_diffusivity has been called:'
end if

apply_mass_diffusivity = .true.

mass_diffusivity_type  = mass_diffusivity_type_arg
C_DDMD                 = C_DDMD_arg

! Check that mass_diffusivity_type is one of the allowable types:
if (my_node .eq. 0) then
   if (mass_diffusivity_type .eq. Yoshizawa_mass_diffusivity) then
      print *, ' mass diffusivity type = Yoshizawa'
      write(lun_info, *) ' mass diffusivity type = Yoshizawa'
   else if (mass_diffusivity_type .eq. dilatation_dependent_mass_diffusivity) then
      print *, ' mass diffusivity type = Dilatation dependent'
      write(lun_info, *) ' mass diffusivity type = Dilatation dependent'
      write(lun_info, *) ' The coefficient C_DDMD = ', C_DDMD
   else
      print *, ' Unrecognized mass diffusivity type = ', mass_diffusivity_type
      stop
   end if
end if

allocate(grad_rho_z  (sr:er, sphi:ephi, nz))
allocate(grad_rho_r  (sr:er, sphi:ephi, nz))
allocate(grad_rho_phi(sr:er, sphi:ephi, nz))
n_words_allocated = n_words_allocated + 3*mr*mphi*nz

! Mass fluxes and mass term:
allocate(mass_fluxes(sr:er, sphi:ephi, nz, 3))
allocate(mass_diff_term  (sr:er, sphi:ephi, nz))
n_words_allocated = n_words_allocated + 2*mr*mphi*nz

if (mass_diffusivity_type .eq. dilatation_dependent_mass_diffusivity) then
   if (.not. allocated(dil)) then
      allocate(dil          (sr:er, sphi:ephi, nz))
      n_words_allocated = n_words_allocated + mr*mphi*nz
   end if
end if

! This is to make sure we have this info. in case the code crashes:
close(lun_info)
open(unit = lun_info, file = 'info.txt', form = 'formatted', status = 'unknown', access = 'append')

end subroutine activate_mass_diffusivity

!----------------------------------------------------------------------------------85   

subroutine add_mass_diffusivity_to_qdot(t, q, qdot, lambda_max)

#ifdef mpi_code
   use mpi
#endif

use grid
use partition_data
use mass_diffusivity
use mass_diffusivity_types
use dof_indices
implicit none
real(8) :: t
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q, qdot
real(8) :: lambda_max ! eigenvalue for time step determination.

! Debug
integer :: ier

! Strain tensor and grad rho calculated internally if needed.
if (mass_diffusivity_type .eq. Yoshizawa_mass_diffusivity) then
#ifdef debug_print
   if (my_node .eq. 0) print *, ' routine add_mass_diff..., node 0: Calling Yoshizawa_mass_flux'
#endif   
   call Yoshizawa_mass_flux           (q, lambda_max) ! subgrid
else if (mass_diffusivity_type .eq. dilatation_dependent_mass_diffusivity) then
#ifdef debug_print
   if (my_node .eq. 0) print *, ' routine add_mass_diff..., node 0: calling dilatation_dependent_mass_flux'
#endif   
   call dilatation_dependent_mass_flux(q, lambda_max) ! subgrid
end if

! Take divergence of the mass fluxes.
#ifdef debug_print
   if (my_node .eq. 0) print *, ' routine add_mass_diff..., node 0: calling mass_diffusivity_term'
#endif
call mass_diffusivity_term

#ifdef debug_print
   if (my_node .eq. 0) print *, ' routine add_mass_diff..., node 0: Returned from mass_diffusivity_term'
#endif

qdot(:,:,:,irho) = qdot(:,:,:,irho) + mass_diff_term(:,:,:)

#ifdef debug_print
   if (my_node .eq. 0) print *, ' routine add_mass_diff..., node 0: Finished adding mass_diff to qdot'
#endif

!call tecplot_mass_diff_term_in_meridional_plane(1, t)

#ifdef mpi_code
   ! This is to make sure that all processors finish what they are doing before
   ! we wrap up.
   ! call mpi_barrier(mpi_comm_world, ier)
   ! call mpi_finalize(ier)
#endif

end subroutine add_mass_diffusivity_to_qdot

!----------------------------------------------------------------------------------85   

subroutine Yoshizawa_mass_flux(q, lambda_max)

! Calculate the Yoshizawa subgrid mass fluxes.

use viscous
use mass_diffusivity
use partition_data
use grid
use dof_indices
use thermal_parameters
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: lambda_max ! eigenvalue for time step determination.

! Locals:
real(8) :: rho, ur, uz, uphi
real(8) :: square, D_u, nu_e
real(8) :: kappa_rho, kappa_rho_1, kappa_rho_2
integer :: iz, iphi, ir
!real(8), parameter :: C_ru1 = 0.22d0
real(8), parameter :: C_ru1 = 0.50d0

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: Yoshizawa_mass_flux: 1st executable'
#endif

if (.not. have_strain_tensor) then
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: Yoshizawa_mass_flux: Calling strain_tensor'
#endif   
   call strain_tensor(q)
end if

if (.not. have_grad_rho) then
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: Yoshizawa_mass_flux: Calling grad_rho'
#endif   
   call grad_rho     (q)
end if 

lambda_max = 0.0d0
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         ! Sij^2:
         square = 2.d0*T_zr(ir,iphi,iz)**2 + 2.d0*T_pz(ir,iphi,iz)**2 + 2.d0*T_pr(ir,iphi,iz)**2 + &
                       T_zz(ir,iphi,iz)**2 +      T_rr(ir,iphi,iz)**2 +      T_pp(ir,iphi,iz)**2
         D_u    = SQRT(2.0d0 * square) ! eq. 127a in Yoshizawa

         kappa_rho_1 = C_ru1**2 * l_grid_squared(ir,iz) * D_u
         ! kappa_rho_2 = C_ru2*ci_squared_initial(ir, iz)/D_u
         ! kappa_rho = kappa_rho_1 + kappa_rho_2
         kappa_rho = kappa_rho_1 ! Will have to resolve later.
         lambda_max = MAX(lambda_max, kappa_rho/min_grid_size(ir,iz)**2)

         mass_fluxes(ir,iphi,iz,r_comp) = kappa_rho * grad_rho_r(ir,iphi,iz)
         mass_fluxes(ir,iphi,iz,z_comp) = kappa_rho * grad_rho_z(ir,iphi,iz)
         mass_fluxes(ir,iphi,iz,p_comp) = kappa_rho * grad_rho_phi(ir,iphi,iz)
      end do
   end do
end do

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: Yoshizawa_mass_flux: Returning'
#endif

end subroutine Yoshizawa_mass_flux

!----------------------------------------------------------------------------------85

subroutine dilatation_dependent_mass_flux(q, lambda_max)

! This is subgrid mass flux.

use mass_diffusivity
use viscous, only: have_strain_tensor
use partition_data
use grid
use dof_indices
use thermal_parameters
use artificial_pressure_module, only: have_dilatation, dil
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: lambda_max ! eigenvalue for time step determination.

! Locals:
real(8) :: rho, ur, uz, uphi
real(8) :: square, D_u, nu_e
real(8) :: kappa_rho, kappa_rho_1, kappa_rho_2
integer :: iz, iphi, ir

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: dilatation_dependent_mass_flux: 1st executable'
#endif

if (.not. have_dilatation) call dilatation(q) ! will use the strain-tensor if we have it.
if (.not. have_grad_rho  ) call grad_rho  (q)

lambda_max = 0.0d0
do iz = 1, nz
   do iphi = sphi, ephi
      do ir = sr, er
         ! Note: Just like for the artificial pressure, we keep only the negative dilatation.
         kappa_rho = -C_DDMD * l_grid_squared(ir,iz) * min(dil(ir,iphi,iz), 0.0d0)         
         lambda_max = MAX(lambda_max, kappa_rho/min_grid_size(ir,iz)**2)

         mass_fluxes(ir,iphi,iz,r_comp) = kappa_rho * grad_rho_r(ir,iphi,iz)
         mass_fluxes(ir,iphi,iz,z_comp) = kappa_rho * grad_rho_z(ir,iphi,iz)
         mass_fluxes(ir,iphi,iz,p_comp) = kappa_rho * grad_rho_phi(ir,iphi,iz)
      end do
   end do
end do

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: dilatation_dependent_mass_flux: Returning'
#endif

end subroutine dilatation_dependent_mass_flux

!----------------------------------------------------------------------------------85

subroutine mass_diffusivity_term

! Computes the mass diffusion term by taking the divergence of the mass fluxes.
! Looked over coding on March 14, 2019.  See notes of that date under the LES
! model tab.

! See notes of Sept. 6, 2018

use mass_diffusivity
use partition_data
use grid
use fargo_or_plotting_shift
use dof_indices

implicit none

! Local:
! These are dimensioned in z-space but will be used in all spaces (since pencils in
! all directions have the same size).
! to do: Make global
real(8), dimension(sr:er, sphi:ephi, nz) :: transposed, deriv, deriv_in_z_space

integer :: ir, iz, iphi

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: mass_diffusivity_term, first executable'
#endif

! The r-component needs to be multiplied by r.
do ir = sr, er
   mass_fluxes(ir,:,:,r_comp) = mass_fluxes(ir,:,:,r_comp) * rgrid(ir)
end do

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: mass_diffusivity_term: Finished mult by r'
#endif

mass_diff_term = 0.0d0

! --------------
! z derivatives:
! --------------
if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
   call pade_diff_bundle(nbundle_z, nz, Ji_z, mass_fluxes(sr,sphi,1,z_comp), deriv_in_z_space)
   mass_diff_term = mass_diff_term + deriv_in_z_space
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: mass_diffusivity_term: Finished z derivatives'
#endif
end if

! ------------------------------
! r derivatives and Fargo terms:
! ------------------------------
if (nr .ne. 1) then
   call transpose_z_to_r (1, mass_fluxes(sr,sphi,1,r_comp), transposed)
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, deriv_in_z_space)
   do ir = sr, er
      mass_diff_term(ir,:,:) = mass_diff_term(ir,:,:) + deriv_in_z_space(ir,:,:)/rgrid(ir)
   end do

   ! Fargo terms:
   if ((add_fargo_extra_operator_now) .and. (nphi .ne. 1)) then
      call transpose_z_to_phi(1, mass_fluxes(sr,sphi,1,r_comp), transposed)
      call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
      call transpose_phi_to_z(1, deriv, deriv_in_z_space)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: mass_diffusivity_term: about to divide by r in r derivatives'
#endif      
      do ir = sr, er
         mass_diff_term(ir,:,:) = mass_diff_term(ir,:,:)+fargo_factor(ir)*deriv_in_z_space(ir,:,:)/rgrid(ir)
      end do
   end if
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: mass_diffusivity_term: Finished r derivatives'
#endif
   
end if

! ----------------
! phi derivatives:
! ----------------
if (nphi .ne. 1) then
   call transpose_z_to_phi(1, mass_fluxes(sr,sphi,1,phi_comp), transposed)
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z(1, deriv, deriv_in_z_space)
#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: mass_diffusivity_term: about to divide by r in phi derivatives'
#endif   
   do ir = sr, er
      mass_diff_term(ir,:,:) = mass_diff_term(ir,:,:) + deriv_in_z_space(ir,:,:)/rgrid(ir)
   end do
end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node 0: mass_diffusivity_term: Returning'
#endif

end subroutine mass_diffusivity_term

!----------------------------------------------------------------------------------85

subroutine grad_rho(q)

! Computes the gradient of the density which can be use for mass diffusivity.
! Looked over coding on March 14, 2019.

use mass_diffusivity
use grid
use partition_data
use fargo_or_plotting_shift
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q

! Local:
! To do: make these global.
real(8), dimension(sr:er, sphi:ephi, nz) :: transposed, deriv
integer :: ir

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: routine grad_rho: first executable'
#endif

if ((nz .ne. 1) .and. (.not. suppress_z_derivatives_when_nz_not_1)) then
   call pade_diff_bundle(nbundle_z, nz, Ji_z, q(sr,sphi,1,1), grad_rho_z)
end if

if (nr .ne. 1) then
   call transpose_z_to_r(1, q(sr,sphi,1,1), transposed)   ! transpose 13
   call pade_diff_bundle(nbundle_r, nr, Ji_r, transposed, deriv)
   call transpose_r_to_z (1, deriv, grad_rho_r)           ! transpose 14
end if

if (nphi .ne. 1) then
   call transpose_z_to_phi(1, q(sr,sphi,1,1), transposed) ! transpose 15
   call pade_diff_periodic(nbundle_phi, nphi, dphi, transposed, deriv)
   call transpose_phi_to_z(1, deriv, grad_rho_phi)    ! transpose 16

   ! Fargo for the r-derivative.  grad_rho_phi is d rho/dphi at the moment:
   if (add_fargo_extra_operator_now) then
      do ir = sr, er
         grad_rho_r(ir,:,:) = grad_rho_r(ir,:,:) + fargo_factor(ir)*grad_rho_phi(ir,:,:)
      end do
   end if

   ! Put in the r factor:
   do ir = sr, er
      grad_rho_phi(ir,:,:) = grad_rho_phi(ir,:,:)/rgrid(ir)
   end do
end if

have_grad_rho = .true.

#ifdef debug_print
if (my_node .eq. 0) print *, ' node 0: routine grad_rho: returning'
#endif

end subroutine grad_rho

!----------------------------------------------------------------------------------85

subroutine tecplot_mass_diff_term_in_meridional_plane(iphi, t)

! This subroutine works in either serial or parallel mode.  The direct access
! and specification of the record number for each data line allows each processor
! to write to the proper record.

use dof_indices
use logical_units
use grid
use thermal_parameters
use partition_data
use basic_state
use mass_diffusivity
use physical_constants
implicit none
integer, intent(in) :: iphi
real(8), intent(in) :: t

! Local:
character(80) :: filename
integer :: ir, iz, irec

write (filename, 4) int(t), t - int(t)
4 format ('mass_diff_meridional_t_', i6.6, f0.4, '.tec')
! recl is the maximum record length.
open (unit = lun_tecplot, file = filename, form = 'formatted', status = 'unknown', access = 'direct', recl = 13*3 + 2)

if (my_node .eq. 0) then
   ! char(10) = new line
   write (lun_tecplot, 1, rec = 1) t, char(13), char(10)
   1  format ('TITLE = "t = ', e12.5, '"', a1, a1)
   write (lun_tecplot, 2, rec = 2) char(13), char(10)
   2  format('VARIABLES = "r", "z", "mass diff"', a1, a1)
   write (lun_tecplot, 3, rec = 3) nr, nz, char(13), char(10)
   3 format('ZONE I=', i4, ',', ' J=',i4, ' DATAPACKING=POINT', a1, a1)
end if

do iz = 1, nz
   do ir = sr, er
      irec = (iz - 1)*nr + ir + 3

      write(lun_tecplot, "(3(1x, e12.5), a1, a1)", rec = irec) rgrid(ir)/AU, zgrid(iz)/AU, &
           mass_diff_term(ir,iphi,iz), char(13), char(10)
   end do
end do
close(lun_tecplot)

if (my_node .eq. 0) print *, ' node 0: Wrote tecplot file ', filename

end subroutine tecplot_mass_diff_term_in_meridional_plane

!----------------------------------------------------------------------------------85
