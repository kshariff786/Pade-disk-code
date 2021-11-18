!----------------------------------------------------------------------------------85

subroutine transpose_z_to_x (nvar, q, q_x_space)

use grid
use partition_data
use transpose_buffer
implicit none
integer :: nvar
real(8), dimension(sx:ex,     sy:ey, nz,   nvar), intent(in)  :: q
real(8), dimension(sy:ey, sz_x:ez_x, nvar, nx  ), intent(out) :: q_x_space

! Local:
integer ivar, iz, iy, ix

#ifdef mpi_code
   call alans_z_to_x (nvar, trans_buffer, q, q_x_space)
#else
   do ivar = 1, nvar
      do iz = 1, nz
         do iy = 1, ny
            do ix = 1, nx
               q_x_space(iy, iz, ivar, ix) = q(ix, iy, iz, ivar)
            end do
         end do
      end do
   end do
#endif

#ifdef debug_print
   if (my_node .eq. 0) print *, ' About to return from transpose_z_to_x'
#endif

end subroutine transpose_z_to_x

!----------------------------------------------------------------------------------85

subroutine transpose_x_to_z (nvar, q_x_space, q)

use grid
use partition_data
use transpose_buffer
implicit none
integer :: nvar
real(8), dimension(sy:ey, sz_x:ez_x, nvar, nx  ), intent(in ) :: q_x_space
real(8), dimension(sx:ex,     sy:ey, nz,   nvar), intent(out) :: q

! Local:
integer ivar, iz, iy, ix

#ifdef debug_print
   if (my_node .eq. 0) print *, ' in transpose_x_to_z.  nvar = ', nvar
#endif

#ifdef mpi_code
   call alans_x_to_z(nvar, trans_buffer, q_x_space, q)
#else
   do ivar = 1, nvar
      do iz = 1, nz
         do iy = 1, ny
            do ix = 1, nx
               q(ix, iy, iz, ivar) = q_x_space(iy, iz, ivar, ix) 
            end do
         end do
      end do
   end do
#endif

#ifdef debug_print
   if (my_node .eq. 0) print *, ' node = 0: About to return from transpose_x_to_z'
#endif

end subroutine transpose_x_to_z
   
!----------------------------------------------------------------------------------85

subroutine transpose_z_to_y (nvar, q, q_y_space)

use grid
use partition_data
use transpose_buffer
implicit none
integer :: nvar
real(8), dimension(sx:ex, sy:ey,     nz,   nvar), intent(in)  :: q
real(8), dimension(sx:ex, sz_y:ez_y, nvar, ny), intent(out) :: q_y_space

! Local:
integer ivar, iz, iy, ix

#ifdef mpi_code
   call alans_z_to_y (nvar, trans_buffer, q, q_y_space)
#else
   do iy = 1, ny
      do ivar = 1, nvar
         do iz = 1, nz
            do ix = 1, nx
               q_y_space(ix, iz, ivar, iy) = q(ix, iy, iz, ivar)
            end do
         end do
      end do
   end do
#endif

#ifdef debug_print
   if (my_node .eq. 0) print *, ' About to return from transpose_z_to_y'
#endif

end subroutine transpose_z_to_y

!----------------------------------------------------------------------------------85

subroutine transpose_y_to_z (nvar, q_y_space, q)

use grid
use partition_data
use transpose_buffer
implicit none
integer :: nvar
real(8), dimension(sx:ex, sz_y:ez_y, nvar, ny  ), intent(in)  :: q_y_space
real(8), dimension(sx:ex, sy:ey,     nz,   nvar), intent(out) :: q

! Local:
integer ivar, iz, iy, ix

#ifdef mpi_code
   call alans_y_to_z (nvar, trans_buffer, q_y_space, q)
#else
   do ivar = 1, nvar
      do iz = 1, nz
         do iy = 1, ny
            do ix = 1, nx
               q(ix, iy, iz, ivar) = q_y_space(ix, iz, ivar, iy)
            end do
         end do
      end do
   end do
#endif

#ifdef debug_print
   if (my_node .eq. 0) print *, ' About to return from transpose_y_to_z'
#endif

end subroutine transpose_y_to_z

!----------------------------------------------------------------------------------85

#ifdef mpi_code

Subroutine alans_z_to_x (nvar, b, z_data, x_data)
  Use MPI
  Use partition_data_for_alan
  Implicit none
  Integer, intent(in) :: nvar
  Real(8), intent(out) :: b(my,mz_x,nvar,mx,0:ng1-1)        ! transpose buffer
  Real(8), intent(in)  :: z_data(mx,my,nz_for_alan,nvar)
  Real(8), intent(out) :: x_data(my,mz_x,nvar,nx)
  Integer :: i,j,k,n,q, ierror

  ! Pack z_data into the buffer b for transmission
  Do q = 0,ng1-1
     Do n = 1,nvar
        Do k = 1,mz_x
           Do j = 1,my                                            ! send all owned y's  (same for sender and receiver)
              Do i = 1,mx                                         ! send all owned x's
                 b(j,k,n,i,q) = z_data(i,j,k+mz_x*q,n)            ! send only receiver's z's
              End Do
           End Do
        End Do
     End Do
  End Do

  Call MPI_ALLTOALL(b, mx*my*mz_x*nvar, MPI_DOUBLE_PRECISION, x_data, &
       mx*my*mz_x*nvar, MPI_DOUBLE_PRECISION, MPI_COMM_xz, ierror)

End Subroutine alans_z_to_x

!----------------------------------------------------------------------------------85  

Subroutine alans_x_to_z(nvar, b, x_data, z_data)
  Use mpi
  Use partition_data_for_alan
  Implicit none
  Integer, intent(in) :: nvar
  Real(8), intent(out) :: b(my,mz_x,nvar,mx,0:ng1-1)        ! transpose buffer
  Real(8), intent(in)  :: x_data(my,mz_x,nvar,nx)
  Real(8), intent(out) :: z_data(mx,my,nz_for_alan,nvar)
  Integer :: i,j,k,n,q, ierror

  Call MPI_ALLTOALL(x_data, mx*my*mz_x*nvar, MPI_DOUBLE_PRECISION, b, mx*my*mz_x*nvar,&
       MPI_DOUBLE_PRECISION, MPI_COMM_xz, ierror)

  ! Unpack the buffer b into z_data
  Do q = 0,ng1-1
     Do n = 1,nvar
        Do j = 1,my                                 ! receive all owned y's  (same for sender and receiver)
           Do k = 1,mz_x
              Do i = 1,mx                               ! receive all owned x's
                 z_data(i,j,k+mz_x*q,n) = b(j,k,n,i,q)  ! receive only sender's z's
              End Do
           End Do
        End Do
     End Do
  End Do

End Subroutine alans_x_to_z

!----------------------------------------------------------------------------------85

Subroutine alans_z_to_y (nvar, b, z_data, y_data)
  Use MPI
  Use partition_data_for_alan
  Implicit none
  Integer, intent(in) :: nvar
  Real(8), intent(out) :: b(mx,mz_y,nvar,my,0:ng2-1)        ! transpose buffer
  Real(8), intent(in)  :: z_data(mx,my,nz_for_alan,nvar)
  Real(8), intent(out) :: y_data(mx,mz_y,nvar,ny)
  Integer :: i,j,k,n,q, ierror

  ! Pack z_data into the buffer b for transmission
  Do q = 0,ng2-1
     Do n = 1,nvar
        Do k = 1,mz_y
           Do j = 1,my                                 ! send all owned y's
              Do i = 1,mx                              ! send all owned x's  (same for sender and receiver)
                 b(i,k,n,j,q) = z_data(i,j,k+mz_y*q,n) ! send only receiver's z's
              End Do
           End Do
        End Do
     End Do
  End Do

  Call MPI_ALLTOALL(b, mx*my*mz_y*nvar, MPI_DOUBLE_PRECISION, y_data, &
       mx*my*mz_y*nvar, MPI_DOUBLE_PRECISION, MPI_COMM_yz, ierror)

End Subroutine alans_z_to_y

!----------------------------------------------------------------------------------85

Subroutine alans_y_to_z (nvar, b, y_data, z_data)
  Use MPI
  Use partition_data_for_alan
  Implicit none
  Integer, intent(in) :: nvar
  Real(8), intent(out) :: b(mx,mz_y,nvar,my,0:ng2-1)        ! transpose buffer
  Real(8), intent(in)  :: y_data(mx,mz_y,nvar,ny)
  Real(8), intent(out) :: z_data(mx,my,nz_for_alan,nvar)
  Integer :: i,j,k,n,q, ierror

  Call MPI_ALLTOALL(y_data, mx*my*mz_y*nvar, MPI_DOUBLE_PRECISION, b, &
       mx*my*mz_y*nvar, MPI_DOUBLE_PRECISION, MPI_COMM_yz, ierror)

  ! Unpack the buffer b into z_data
  Do q = 0,ng2-1
     Do n = 1,nvar
        Do k = 1,mz_y
           Do j = 1,my                                            ! receive all owned y's
              Do i = 1,mx                                         ! receive all owned x's  (same for sender and receiver)
                 z_data(i,j,k+mz_y*q,n) = b(i,k,n,j,q)            ! receive only sender's z's
              End Do
           End Do
        End Do
     End Do
  End Do

End Subroutine alans_y_to_z

!----------------------------------------------------------------------------------85  

subroutine make_communicators_for_alans_transposes

! Form the MPI communicators for the x- and y-transposes, namely,
! mpi_comm_xz and mpi_comm_yz.

! In addition we determine the grid point ranges assigned to each processor in each
! of the phases.  These are denoted sx, er, etc.

use mpi
use partition_data
use partition_data_for_alan, only: mx, mz_x, my, mz_y, ng1, ng2, mpi_comm_xz, mpi_comm_yz
use boundary_conditions
implicit none
integer, allocatable, dimension(:) :: group_ranks
integer :: px, py
logical :: found
integer :: mpi_group_world, mpi_group_scratch, mpi_comm_scratch, ierror

integer :: sx, sz_x, sy, sz_y

! Code copied from Alan's stellar box code.

allocate(group_ranks(0 : max(ng1, ng2)-1))

! make a group out of the full set of processors.  This is needed below.
Call mpi_comm_group(mpi_comm_world, mpi_group_world, ierror)

! These are not needed and were put only so that the compiler did not
! issue a warning about a possible unitialization.
sx   = 1
sz_x = 1

Do py = 0, ng2 - 1
   found = .false.
   Do px = 0, ng1 - 1
      group_ranks(px) = px + ng1*py
      If (group_ranks(px) == my_node) then
         sx   = 1 + mx*px      ! first x-node when in y- or z-space
         sz_x = 1 + mz_x*px    ! first z-node when in x-space
         ! stz_x = 1 + mtz_x*px  ! first z-node when in ray-tracing x-space
         found = .true.
      End If
   End Do
   call mpi_group_incl (mpi_group_world, ng1, group_ranks, mpi_group_scratch, ierror)
   call mpi_comm_create(mpi_comm_world, mpi_group_scratch, mpi_comm_scratch, ierror)
   If (found) mpi_comm_xz = mpi_comm_scratch
End Do

! These are not needed and were put only so that the compiler did not
! issue a warning about a possible unitialization.
sy   = 1
sz_y = 1

Do px = 0,ng1-1
   found = .false.
   Do py = 0,ng2-1
      group_ranks(py) = px + ng1*py
      If (group_ranks(py)==my_node) then
         sy   = 1 + my*py      ! first y-node when in x- or z-space
         sz_y = 1 + mz_y*py    ! first z-node when in y-space
         ! stz_y = 1 + mtz_y*py  ! first z-node when in ray_tracing y-space
         found = .true.
      End If
   End Do
   call mpi_group_incl(mpi_group_world, ng2, group_ranks, mpi_group_scratch, ierror) 
   call mpi_comm_create(mpi_comm_world, mpi_group_scratch, mpi_comm_scratch, ierror)
   If (found) mpi_comm_yz = mpi_comm_scratch
End Do

deallocate(group_ranks)

ex   = sx     + (mx     - 1)
ez_x = sz_x   + (mz_x   - 1)
ey   = sy   + (my   - 1)
ez_y = sz_y + (mz_y - 1)

end subroutine make_communicators_for_alans_transposes

#endif

!----------------------------------------------------------------------------------85
