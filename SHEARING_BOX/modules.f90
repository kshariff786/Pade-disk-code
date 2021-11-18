!----------------------------------------------------------------------------------85

! Use these to improve code readability.
module dof_indices
   implicit none
   integer, parameter :: irho = 1, xmom = 2, zmom = 3, ymom = 4, ener = 5
end module dof_indices

module grid
   implicit none
   logical :: gridded, periodic_z
   integer, parameter :: ndof = 5
   integer :: nx, nz, ny

   ! Grid sizes.  Only z direction can have a non-uniform mesh.
   real(8) :: dx, dy, dx_inv, dy_inv
   real(8), dimension(:), allocatable :: dz
   real(8) :: dz_periodic
   
   ! Grid point locations:
   real(8), dimension(:), allocatable :: xgrid, zgrid, ygrid

   ! Inverse Jacobian for z direction:
   real(8), dimension(:), allocatable :: Ji_z

   ! Domain.
   real(8) :: xmin, xmax, zmin, zmax, ymin, ymax

   ! Weights for trapezoidal rule integration:
   real(8), dimension(:), allocatable :: trap_weight_z, trap_weight_y

   ! This flag is needed to run an nz = 1 case, i.e., (x, y) plane case in
   ! parallel.  We fake it by using nz = 2 and suppress z derivatives and smoothing.
   logical :: suppress_z_derivatives_when_nz_not_1
   integer :: nz_actual
end module grid 

module logical_units
   implicit none
   integer, parameter :: lun_general_purpose = 1, lun_tecplot = 2
end module logical_units

module math_constants
   implicit none
   real(8), parameter :: athird = 1.d0/3.d0, asixth = 1.d0/6.d0
end module math_constants

module pade_filter
   logical :: apply_pade_filter
   ! Three characters ('tau' or 'eps').
   ! 'eps' : use the specified eps_filter
   ! 'tau' : determine eps_filter from tau.
   character(3) :: eps_or_tau
   real(8) eps_filter, tau_filter
end module pade_filter   

module partition_data
   logical :: partitioned
   integer :: num_nodes, my_node
   ! Range of indices held by processors.
   integer :: sx, sz_x, sy, sz_y
   integer :: er, ez_x, ey, ez_y
   integer :: mr, mz_x, my, mz_y
   integer :: nbundle_z, nbundle_x, nbundle_y
   integer :: mpi_comm_xz, mpi_comm_yz
   integer :: nz_for_alan ! Since I don't want to use grid in Alan's routines.
   integer :: ng1, ng2   
end module partition_data

module physical_constants
   implicit none
   real(8), parameter :: Gconst = 6.673d-8 ! cm^3 g^-1 s^-2

   ! Stefan-Boltzmann:
   real(8), parameter :: sigma_SB = 5.6705d-5 ! erg s^-1 cm^-2 K^-4
   ! Boltzmann:
   real(8), parameter :: kB       = 1.3807d-16 ! erg K^-1
   real(8), parameter :: Msolar   = 1.989d33 ! gm
   real(8), parameter :: Rsolar   = 69.57d9  ! cm
   real(8), parameter :: AU       = 1.496d13 ! cm
   real(8), parameter :: parsec   = 3.0857d18 ! cm
   real(8), parameter :: kilometer_per_sec = 1000.d0*100.d0 ! cm/sec
   real(8), parameter :: Rgas     = 3.5871d7  ! cm^2 s^-2 K^-1
   real(8), parameter :: year     = 3.16d7  ! s
   ! Mass of atomic hydrogen:
   real(8), parameter :: mH       = 1.6735d-24 ! gm
   real(8), parameter :: mH2      = 2.d0 * mH  ! gm
end module physical_constants

! Flow field array.
module q_array
   implicit none
   real(8), allocatable, dimension(:, :, :, :) :: q
end module q_array

! Similarly, this module was created to avoid dynamic storage allocations and
! deallocations.
module transposes_of_q_and_qdot
   real(8), allocatable, dimension(:, :, :, :) :: q_x_space, qdot_x_space
   real(8), allocatable, dimension(:, :, :, :) :: q_y_space, qdot_y_space
end module transposes_of_q_and_qdot

module transpose_buffer
   ! allocate(trans_buffer(sx:sx+mx-1, sy:sy+my-1, nz, nvar))

   ! Declaration in z_to_x and x_to_z:
   ! b(my,mz_x,nvar,mx,0:ng1-1)  mz_x*ng1 = nz

   ! Declaration in z_to_y and y_to_z
   ! b(mx,mz_y,nvar,my,0:ng1-1)  mz_y*ng2 = nz

   real(8), allocatable, dimension(:, :, :, :) :: trans_buffer
end module transpose_buffer
   
! This module was created to avoid dynamic storage allocations and deallocations
! everytime rk4 is called.  These arrays are allocated in subroutine initialize.
module rk4_arrays
   ! Accumulator, next argument to rhs, rhs array.
   real(8), allocatable, dimension(:, :, :, :) :: q_accum, q_next_arg, qdot
end module rk4_arrays
   
module thermal_parameters
   real(8) :: gamma, gm1
   logical :: isothermal
   ! Initial isothermal speed of sound squared.  Dimensions are (ny, nz).
   ! Isothermality is imposed by setting p = rho * ci_squared_initial(r, z)
   real(8), allocatable, dimension(:, :) :: ci_squared_initial
end module thermal_parameters

!----------------------------------------------------------------------------------85
