!----------------------------------------------------------------------------------85

subroutine apply_plotting_shifts(q, t)

! This routine is used for runs in which a portion of the phi (azimuthal) direction
! is simulated with periodic boundary conditions and you want the plotting to rotate
! with the Keplerian flow, say.  
  

use grid
use partition_data
use fargo_trick_or_plotting_shift
implicit none
real(8), dimension(sr:er, sphi:ephi, nz, ndof) :: q
real(8) :: t

! Local:
integer :: ir, k
complex(8) :: ci = CMPLX(0.0, 1.0d0)
real(8) :: desired_phi_shift

desired_phi_shift = Omega0_plotting_shift * (t - t_previous_plotting_shift)

do ir = sr, er
   do k = 1, nphi/2+1
      complex_shift_factor(k, ir) = EXP(-ci * (k-1) * desired_phi_shift)
   end do
end do

#ifdef fftw
   ! call apply_real_shifts_using_fftw(q)
   call apply_real_shifts_using_fftw_old(q) ! this can actually be cheaper.
#else
   call apply_real_shifts_using_rogallo_fft(q)
#endif

t_previous_plotting_shift = t

end subroutine apply_plotting_shifts

!----------------------------------------------------------------------------------85
