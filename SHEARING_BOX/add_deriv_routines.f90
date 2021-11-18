!----------------------------------------------------------------------------------85

subroutine add_z_derivatives(q, qdot)

! This uses the straight Pade scheme for differentiation.

use grid
use partition_data
use dof_indices, only: irho, ymom, zmom, xmom, ener
use thermal_parameters
use boundary_condition_types
use boundary_conditions
use gravity
use artificial_pressure_module

implicit none
real(8), dimension(sx:ex, sy:ey, nz, ndof) :: q, qdot

! Locals:
! Later see if these can be made static.
real(8), dimension(sx:ex, sy:ey, nz) :: F, Fd
real(8), dimension(sx:ex, sy:ey, nz) :: uz, pressure, d_uz_dz

integer :: iy, iz, ix
real(8) :: rho_2_q2, rho_q2, gravity_term
real(8) :: c_sound, uz1, lambda_max_local, gamma_sound_speed, u_plus_c

! To store boundary flux derivatives for boundary condition modification.
real(8), dimension(sx:ex, sy:ey, ndof) :: dF_bottom, dF_top

! The order of calculation is necessitated by the indexing order of "q" which
! in turn is dictated by Alan's transpose routines.

#ifdef debug_print
   if (my_node .eq. 0) print *, ' first executable in add_z_derivatives'
#endif

if (isothermal) then
   gamma_sound_speed = 1.0d0
else
   gamma_sound_speed = gamma
end if

! Store uz which is needed below:
uz(:, :, :) = q(:, :, :, zmom) / q(:, :, :, irho)

! Put pressure into an array which is needed below:
if (isothermal) then
   do iz = 1, nz
      do iy = sy, ey
         do ir = sx, ex
            pressure(ix, iy, iz) = q(ix, iy, iz, irho) * ci_squared_initial(ix, iz)
         end do
      end do
   end do
else
   do iz = 1, nz
      do iy = sy, ey
         do ir = sx, ex
            pressure(ix, iy, iz) = q(ix, iy, iz, ener) * gm1 ! From internal energy.
         end do
      end do
   end do
end if

if (apply_artificial_pressure) pressure = pressure + p_art

! Mass equation:
! ~~~~~~~~~~~~~~
F(:, :, :) = q(:, :, :, zmom) ! rho*uz
call pade_diff_z(nbundle_z, F, Fd)

! Store for later modification by boundary conditions:
dF_bottom(:, :, irho) = Fd(:, :, 1 )
dF_top   (:, :, irho) = Fd(:, :, nz)

qdot(:, :, :, irho) = qdot(:, :, :, irho) - Fd(:, :, :)

! y-momentum equation:
! ~~~~~~~~~~~~~~~~~~~~
F(:, :, :) = q(:, :, :, ymom) * uz(:, :, :)
call pade_diff_z(nbundle_z, F, Fd)

! Store for possible later modification by non-reflective boundary conditions:
dF_bottom(:, :, ymom) = Fd(:, :, 1 )
dF_top   (:, :, ymom) = Fd(:, :, nz)

qdot(:, :, :, ymom) = qdot(:, :, :, ymom) - Fd(:, :, :)

! x-momentum equation:
! ~~~~~~~~~~~~~~~~~~~~~~~~~
F(:, :, :) = q(:, :, :, rmom) * uz(:, :, :)
call pade_diff_z(nbundle_z, F, Fd)
dF_bottom(:, :, rmom) = Fd(:, :, 1 )
dF_top   (:, :, rmom) = Fd(:, :, nz)
qdot(:, :, :, rmom) = qdot(:, :, :, rmom) - Fd(:, :, :)

! Vertical momentum equation:
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
F(:, :, :) = pressure(:,:,:) + q(:,:,:,zmom) * uz(:,:,:)
call pade_diff_z(nbundle_z, F, Fd)
dF_bottom(:, :, zmom) = Fd(:, :, 1 )
dF_top   (:, :, zmom) = Fd(:, :, nz)
qdot(:, :, :, zmom) = qdot(:, :, :, zmom) - Fd(:, :, :)

! Internal energy equation:
! ~~~~~~~~~~~~~~~~~~~~~~~~~
if (.not. isothermal) then
   ! d/dz(e_int uz):
   F(:, :, :) = q(:, :, :, ener)*uz(:, :, :)
   call pade_diff_z(nbundle_z, F, Fd)
   
   ! + p div_u term on LHS:
   call pade_diff_z(nbundle_z, uz, d_uz_dz)

   ! Combine the two terms into Fd:
   Fd(:,:,:) = Fd(:,:,:) + pressure(:,:,:)*d_uz_dz(:,:,:)            

   ! Store for possible later modification by boundary conditions:
   dF_bottom(:, :, ener) = Fd(:, :, 1 )
   dF_top   (:, :, ener) = Fd(:, :, nz)

   ! Add to qdot:
   qdot(:, :, :, ener) = qdot(:, :, :, ener) - Fd(:, :, :)
end if

! Non-reflective boundary conditions at the top and bottom boundaries.
! qdot is modified.
if  (zmin_BC .eq. non_reflective) then
   iz = 1
   if (isothermal) then
      call z_isothermal_non_reflective_bc(iz, q, uz, dF_bottom, qdot)
   else
      call z_adiabatic_non_reflective_bc (iz, q, uz, dF_bottom, qdot)
   end if
end if

if  (zmax_BC .eq. non_reflective) then
   iz = nz
   if (isothermal) then
      call z_isothermal_non_reflective_bc(iz, q, uz, dF_top, qdot)
   else
      call z_adiabatic_non_reflective_bc (iz, q, uz, dF_top, qdot)      
   end if
end if
! Done with non-reflective BC 
   
#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' returning from add_z_derivatives'
   end if   
#endif

return
end subroutine add_z_derivatives

!----------------------------------------------------------------------------------85

subroutine add_x_derivatives (q, qdot)

! Adds x derivative terms to qdot.
! q and qdot should be in z space.
! We will perform a transpose to x space and then back to z space here.

use grid, only: ndof, nx, nz, Ji_x, xgrid, dy, Ji_x, dx
use partition_data
use dof_indices, only: irho, xmom, ymom, zmom, ener
use thermal_parameters
use transposes_of_q_and_qdot
use artificial_pressure_module
use boundary_condition_types
use boundary_conditions
use thermal_parameters
implicit none

real(8), dimension(sx:ex, sy:ey, nz, ndof), intent(in)  :: q
real(8), dimension(sx:ex, sy:ey, nz, ndof), intent(out) :: qdot

! Local:
integer :: ix, iz, iy, idof

! Again, see later if these can be made static.
real(8), dimension(sy:ey, sz_x:ez_x, ndof, nx) :: F, dFdx

! Needed since using the internal energy:
real(8), dimension(sy:ey, sz_x:ez_x, nx) :: r_ur, div_u_term

! This is used so often that I decided to make it an array:
real(8), dimension(sy:ey, sz_x:ez_x, nx) :: ur

real(8) :: c_sound, lambda_max_local, gamma_sound_speed

! Pressure associated locals:
real(8), dimension(sy:ey, sz_x:ez_x, nx):: pressure
real(8) :: rho_2_q2, rho_q2

#ifdef debug_print
   ! Note: my_node = 0 for the serial code so this works both cases.
   if (my_node .eq. 0) then
      print *, ' add_x_derivatives has been called'
   end if
#endif

! Perform transposes:
call transpose_z_to_x (ndof, q,    q_x_space   )
call transpose_z_to_x (ndof, qdot, qdot_x_space)

if (isothermal) then
   gamma_sound_speed = 1.0d0
else
   gamma_sound_speed = gamma
end if

! Store pressure which is needed below:
if (isothermal) then
   do ir = 1, nx
      do iz = sz_x, ez_x
         do iy = sy, ey
            pressure(iy, iz, ir) = q_x_space(iy, iz, irho, ir) * ci_squared_initial(ix, iz)
         end do
      end do
   end do
else
   do ir = 1, nx
      do iz = sz_x, ez_x
         do iy = sy, ey
            pressure(iy, iz, ir) = q_x_space(iy, iz, ener, ir) * gm1 ! from internal energy
         end do
      end do
   end do
end if

! Artificial pressure using dil in r-space:
if (apply_artificial_pressure) then
   call transpose_z_to_x (1, p_art, p_art_x_space)
   do ir = 1, nx
      do iz = sz_x, ez_x
         do iy = sy, ey
            pressure(iy, iz, ir) = pressure(iy, iz, ir) + p_art_x_space(iy, iz, ir)
         end do
      end do
   end do
end if

! Store uy:
do ir = 1, nx
   do iz = sz_x, ez_x
      do iy = sy, ey
          uy(iy, iz, ir) = q_x_space(iy, iz, rmom, ir) / q_x_space(iy, iz, irho, ir)
      end do
   end do
end do

! Compute fluxes.  On Oct 25, 2018 I put in the pressure into the flux in order to
! implement the non-reflective condition.  This means that I need a p/r term on the rhs.
do ir = 1, nx
   do iz = sz_x, ez_x
      do iy = sy, ey
         F(iy,iz,irho,ir) = rgrid(ix)* q_x_space(iy,iz,rmom,ir)
         F(iy,iz,ymom,ir) = rgrid(ix)* q_x_space(iy,iz,ymom,ir)*uy(iy,iz,ir)
         F(iy,iz,rmom,ir) = rgrid(ix)*(q_x_space(iy,iz,rmom,ir)*uy(iy,iz,ir) + pressure(iy,iz,ir))
         F(iy,iz,zmom,ir) = rgrid(ix)* q_x_space(iy,iz,zmom,ir)*uy(iy,iz,ir)
      end do
   end do
end do

if (.not. isothermal) then
   do ir = 1, nx
      do iz = sz_x, ez_x
         do iy = sy, ey
            F(iy, iz, ener,  ir) = q_x_space(iy,iz,ener,ir) * rgrid(ix) * uy(iy,iz,ir)
            r_uy(iy, iz, ir) = uy(iy,iz,ir) * rgrid(ix)
         end do
      end do
   end do
end if

! Why do I have this?
!   do ir = 1, nx
!      do iz = sz_x, ez_x
!         do iy = sy, ey
!            F(iy,iz,ener,ir) = (q_x_space(iy,iz,ener,ir) + pressure(iy,iz,ir)) * rgrid(ix) * ur(iy,iz,ir)
!         end do
!      end do
!   end do
!end if

! Differentiate the fluxes:
! Note that for the isothermal case we are differentiating one more dof
! than needed.  This can't be helped due to the indexing of F we have chosen.   
call pade_diff_bundle(my*mz_x*ndof, nx, Ji_x, F, dFdr)

! Non-reflective boundary conditions at the left and right boundaries.
! dFdr is modified:
if  (rmin_BC .eq. non_reflective) then
   if (isothermal) then
      ir = 1
      ! print *, ' calling r_isothermal_non_reflective_bc'
      call r_isothermal_non_reflective_bc(ix, q_x_space, uy, dFdr)
   end if
end if

if  (rmax_BC .eq. non_reflective) then
   if (isothermal) then
      ir = nx
      call r_isothermal_non_reflective_bc(ix, q_x_space, uy, dFdr)
   end if
end if
! Done with non-reflective BC 
   
   
if (.not.isothermal) then
   call pade_diff_bundle(my*mz_x, nx, Ji_x, r_uy, div_u_term)
end if   

do ir = 1, nx
   do idof = 1, ndof
      do iz = sz_x, ez_x
         do iy = sy, ey
            qdot_x_space(iy, iz, idof, ir) = qdot_x_space(iy, iz, idof, ir) - &
                                       dFdr(iy, iz, idof, ir)/rgrid(ix)
         end do
      end do
   end do
end do

! - p div.u term in internal energy equation:
if (.not.isothermal) then
   do ix = 1, nx
      do iz = sz_x, ez_x
         do iy = sy, ey
            qdot_x_space(iy, iz, ener, ir) = qdot_x_space(iy, iz, ener, ix) - &
                          pressure(iy, iz, ix) * div_u_term(iy, iz, ix)
         end do
      end do
   end do
end if

! Transpose qdot back to z space:
call transpose_x_to_z (ndof, qdot_x_space, qdot)

return
end subroutine add_x_derivatives

!----------------------------------------------------------------------------------85

subroutine add_y_derivatives(q, qdot)

! Adds y derivative terms to qdot.
! Both q and qdot should have z-space indexing.  We will perform a transpose here.

! Note: The extra Fargo terms are treated in subroutine add_fargo_extra_operator_terms.

use transposes_of_q_and_qdot
use grid
use partition_data
use dof_indices, only: irho, zmom, rmom, ymom, ener
use thermal_parameters
use fargo_or_plotting_shift
use artificial_pressure_module
implicit none

real(8), dimension(sx:ex, sz_y:ez_y, ndof, ny), intent(in)    :: q
real(8), dimension(sx:ex, sz_y:ez_y, ndof, ny), intent(inout) :: qdot

! Local:
integer :: ir, iz, iy, idof

! Again, see later if these can be made static.
real(8), dimension(sx:ex, sz_y:ez_y, ndof, ny) :: F, dF
real(8) :: c_sound, lambda_max_local, gamma_sound_speed

! Pressure associated locals:
real(8), dimension(sx:ex, sz_y:ez_y, ny) :: p, dp_dy
real(8) :: rho_2_q2, rho_q2

real(8), dimension(sx:ex, sz_y:ez_y, ny) :: uy
! Residual azimuthal velocity for Fargo trick:
real(8), dimension(sx:ex, sz_y:ez_y, ny) :: uy_prime

! Needed for the internal energy.
real(8), dimension(sx:ex, sz_y:ez_y, ny) :: div_u_term

! For Fargo trick:
!real(8) :: lambda_max_local_without_fargo

! Artificial pressure related local:
real(8) :: beta_ap

! du and du_dy for the shear advection terms:
real(8), dimension(sx:ex, sz_y:ez_y, 3, ny) :: u, du

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' add_y_derivatives has been called'
   end if
#endif

! Transposes:
#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' Calling transpose_z_to_y for q in add_y_derivatives'
   end if
#endif
call transpose_z_to_y (ndof, q, q_y_space)
#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' Calling transpose_z_to_y for qdot in add_y_derivatives'
   end if
#endif   
call transpose_z_to_y (ndof, qdot, qdot_y_space)   

if (isothermal) then
   gamma_sound_speed = 1.0d0
else
   gamma_sound_speed = gamma
end if

#ifdef debug_print
   if (my_node .eq. 0) print *, ' add_y_derivatives: about to compute pressure'
#endif

! Store pressure which is needed below:
if (isothermal) then
   do iy = 1, ny
      do iz = sz_y, ez_y
         do ix = sx, ex
            p(ix, iz, iy) = q_y_space(ix, iz, irho, iy) * ci_squared_initial(ix, iz)
         end do
      end do
   end do
else
   do iy = 1, ny
      do iz = sz_y, ez_y
         do ir = sx, er
            p(ix, iz, iy) = q_y_space(ix, iz, ener, iy) * gm1
         end do
      end do
   end do
end if

! Add artificial pressure to the pressure:
if (apply_artificial_pressure) then
   call transpose_z_to_y (1, p_art, p_art_y_space)
   do iy = 1, ny
      do iz = sz_y, ez_y
         do ir = sx, er
            p(ix, iz, iy) = p(ix, iz, iy) + p_art_y_space(ix, iz, iy)
         end do
      end do
   end do
end if


! Compute fluxes:
do iy = 1, ny
   do iz = sz_y, ez_y
      do ir = sx, er
         F(ix, iz, irho, iy) =  q_y_space(ix, iz, irho, iy) * uy(ix, iz, iy)
         F(ix, iz, ymom, iy) =  q_y_space(ix, iz, ymom, iy) * uy(ix, iz, iy) / rgrid(ix)
         F(ix, iz, rmom, iy) =  q_y_space(ix, iz, rmom, iy) * uy(ix, iz, iy)
         F(ix, iz, zmom, iy) =  q_y_space(ix, iz, zmom, iy) * uy(ix, iz, iy)

          c_sound = SQRT(gamma_sound_speed * p(ix, iz, iy) / q_y_space(ix, iz, irho, iy))
          lambda_max_local = (abs(uy(ix,iz,iy)) + c_sound) / (rgrid(ix)*dy)
          lambda_max = MAX(lambda_max, lambda_max_local)
      end do
   end do
end do

! Energy fluxes:
if (.not. isothermal) then
   do iy = 1, ny
      do iz = sz_y, ez_y
         do ir = sx, ex
            F(ix, iz, ener, iy) = q_y_space(ix, iz, ener, iy) * uy(ix, iz, iy)
         end do
      end do
   end do
else
   do iy = 1, ny
         do iz = sz_y, ez_y
            do ir = sx, er
               F(ix, iz, ener, iy) = (q_y_space(ix,iz,ener,iy) + p(ix,iz,iy)) * uy(ix,iz,iy)
            end do
         end do
      end do
   end if
end if 

! Differentiate the fluxes.  Note: We have an extra un-needed one for the isothermal case:
call pade_diff_periodic(mr*mz_y*ndof, ny, dy, F(sx, sz_y, 1, 1), dF(sx, sz_y, 1, 1))

do iy = 1, ny
   do iz = sz_y, ez_y
      do ir = sx, er
         qdot_y_space(ix,iz,irho,iy) = qdot_y_space(ix,iz,irho,iy) - dF(ix,iz,irho,iy)/rgrid(ix)
         qdot_y_space(ix,iz,rmom,iy) = qdot_y_space(ix,iz,rmom,iy) - dF(ix,iz,rmom,iy)/rgrid(ix)
         qdot_y_space(ix,iz,zmom,iy) = qdot_y_space(ix,iz,zmom,iy) - dF(ix,iz,zmom,iy)/rgrid(ix)
         qdot_y_space(ix,iz,ymom,iy) = qdot_y_space(ix,iz,ymom,iy) - dF(ix,iz,ymom,iy) ! No r
      end do
   end do
end do

if (.not. isothermal) then
   do iy = 1, ny
      do iz = sz_y, ez_y
         do ir = sx, er
            qdot_y_space(ix,iz,ener,iy) = qdot_y_space(ix,iz,ener,iy) - dF(ix,iz,ener,iy)/rgrid(ix)
         end do
      end do
   end do
end if

! - p div u term for internal energy equation:
if (.not.isothermal) then
   call pade_diff_periodic(mr*mz_y, ny, dy, uy, div_u_term)
   do iy = 1, ny
      do iz = sz_y, ez_y
         do ir = sx, er
            qdot_y_space(ix,iz,ener,iy) = qdot_y_space(ix,iz,ener,iy) - &
                 p(ix,iz,iy)/rgrid(ix)*div_u_term(ix,iz,iy)
         end do
      end do
   end do
end if                           

! Y pressure derivative for angular momentum eq.:
call pade_diff_periodic(mr*mz_y, ny, dy, p(sx, sz_y, 1), dp_dy(sx, sz_y, 1))

do iy = 1, ny
   do iz = sz_y, ez_y
      do ir = sx, er
         qdot_y_space(ix,iz,ymom,iy) = qdot_y_space(ix,iz,ymom,iy) - dp_dy(ix, iz, iy)
      end do
   end do
end do

! Advection due to the mean shear:
do iy = 1, ny
   do iz = sz_y, ez_y
      do idir = 1, 3
         do ix = sx, ex
            u(ix, iz, iy) = q_y_space(ix, iz, idir, iy) / q_y_space(ix, iz, irho, iy)
         end do
      end do
   end do
end do

call pade_diff_periodic(mx*mz_y*3, ny, dy, u(sx, sz_y, 1, 1), du(sx, sz_y, 1, 1))

do iy = 1, ny
   do iz = sz_y, ez_y
      do idir = 1, 3
         do ix = sx, ex
            qdot_y_space(ix,iz,xmom,iy) = qdot_y_space(ix,iz,xmom,iy) - Shear*xgrid(ix)*du(ix,izx
            u(ix, iz, iy) = q_y_space(ix, iz, idir, iy) / q_y_space(ix, iz, irho, iy)
         end do
      end do
   end do
end do




end subroutine add_y_derivatives

!----------------------------------------------------------------------------------85


