!**********************************************************************************85
! Disk code in cylindrical coordinates (z, r, phi). 
! Degrees of freedom in the "q" array are:
!        1: Density: rho
!        2: Radial momentum: rho*ur
!        3: Vertical momentum: rho*uz
!        4: Angular momentum: rho*u_phi*r
!        5: Energy per unit volume
!        If internal_energy = .true. then we use e = p/(gamma - 1)
!        Otherwise we use e = p/(gamma - 1) + 1/2 rho q^2.
!        If isothermal = .true. we don't use the energy equation at all and instead
!        use p = rho * ci_squared_initial(ir, iz) where ci_squared_initial is a
!        prescribed isothermal sound speed.

! Note: For the Fargo option, we must use the internal energy as explained in the
!       Ap. J. manuscript.

! For the isothermal option (isothermality = true) we don't need the energy
! degree of freedom but to keep the coding simple, it is carried anyway without
! harm.

! Nov. 16, 2017: Removed the gravitational potential from the definition of energy
!                so now we have gravity work terms in the energy equation.

! Spatial differencing : 4th order Pade.
! Time stepping        : 4th order Runge Kutta.

! For reference the pencil directions are:

! Alan      Me
! x    <--> r
! y    <--> phi
! z    <--> z

! Alan's x-space (y, z, q, x) ---> Our r space   (phi, z, q, r)
! Alan's y-space (x, z, q, y) ---> Our phi space (r, z, q, phi)
! Alan's z-space (x, y, z, q) ---> Our z space   (r, phi, z, q)

! Dimensioning of each pencil:
! q           (sr:er, sphi:ephi, nz, ndof)           : in z space
! q_r_space   (sphi:ephi, sz_r:ez_r, ndof, nr)
! q_phi_space (sr:er, sz_phi:ez_phi, ndof, nphi)

! Things to do:
!    (1) Make ci_squared locally dimensioned to each processor.
!    (2) Make F and Fd global.
!
!----------------------------------------------------------------------------------85

! Main program
program disk

#ifdef mpi_code
   use mpi
#endif
use logical_units, only : lun_general_purpose
use partition_data

implicit none
integer :: i_run_type, ier

#ifdef mpi_code
   call mpi_init (ier)
   call mpi_comm_rank (mpi_comm_world, my_node,   ier)
   call mpi_comm_size (mpi_comm_world, num_nodes, ier)
#else
   print *, ' Serial run'
   my_node = 0 ! For a serial run.
#endif

! This subroutine does the broadcast:
call read_namelist_run_type (lun_general_purpose, i_run_type)

if (i_run_type .eq. 0) then
   call user_application
else if (i_run_type .eq. 1) then
   call euler1d_tests ! in app_euler1d_tests.f90
else if (i_run_type .eq. 2) then
   call hydrostatic_test    ! in app_hydrostatic_test.f90
else if (i_run_type .eq. 3) then
   ! in app_homentropic_solid_body_rotation_test.f90
   call homentropic_solid_body_rotation_test
else if (i_run_type .eq. 4) then
   print *, ' This run type no longer exists'
   call terminate_with_no_save(1)
   ! call axi_vertical_shear_instability
else if (i_run_type .eq. 5) then
   call vertically_integrated_disk
else if (i_run_type .eq. 6) then
   call vsi_3D
else if (i_run_type .eq. 7) then
   print *, ' This run type no longer exists'
   call terminate_with_no_save(1)   
   ! call strong_vortex
else if (i_run_type .eq. 8) then
   print *, ' This run type no longer exists'
   call terminate_with_no_save(1)   
   ! call two_vortices
else if (i_run_type .eq. 9) then
   call vortex_pair
!else if (i_run_type .eq. 10) then
!   call Cassen_Moosman
else if (i_run_type .eq. 11) then
   call taylor_couette
else if (i_run_type .eq. 12) then
   call app_zvi ! located in app_zvi.f90
else
   if (my_node .eq. 0) print *, ' node 0: Unknown i_run_type = ', i_run_type
end if

end program disk

!----------------------------------------------------------------------------------85
