!----------------------------------------------------------------------------------85

! Test of errors in Shu and Osher telescoping sum derivative.  In particular,
! why it converges to second order at the boundaries.

! Number of intervals = # of cell centers.
integer :: nz, ires, nres, nz1

integer :: iz
real(8) :: pi, amp, dz
real(8), allocatable, dimension(:) :: zcen, F, Fd, exact_deriv, cc_nz
real(8), allocatable, dimension(:) :: zfaces, H, Hp, H_exact, Hp_exact, cc_nzp1
real(8), allocatable, dimension

pi = 4.0d0 * atan(1.0d0)

write (6, "(' enter # of resolution levels and nz1--->', $)")
read (5, *) nres, nz1

open (unit = 9,  file = 'err_Shu_Osher_vs_nz.dat', form = 'formatted', status = 'unknown')
open (unit = 10, file = 'err_straight_Pade_vs_nz.dat',   form = 'formatted', status = 'unknown')
open (unit = 11, file = 'err_Hp_vs_nz.dat',   form = 'formatted', status = 'unknown')
open (unit = 12, file = 'err_Hex_p_vs_nz.dat',   form = 'formatted', status = 'unknown')

do ires = 1, nres
   nz = nz1 * 2**(ires - 1)
   print *, ' nz = ', nz
   allocate (zcen(nz), F(nz), Fd(nz), exact_deriv(nz), cc_nz(nz))
   allocate (zfaces(nz+1), H(nz+1), Hp(nz+1), H_exact(nz+1), Hp_exact(nz+1), &
             cc_nzp1(nz+1))


   ! Compute q:
   ! ~~~~~~~~~~
   dz = 2.d0*pi / nz
   do iz = 1, nz
      ! Cell centers:
      zcen(iz)    = 0.5d0*dz + (iz - 1)*dz
      amp = 2.d0 / dz * sin(0.5d0*dz)
      F(iz)           =   amp * cos(zcen(iz))
      exact_deriv(iz) = - amp * sin(zcen(iz))
   end do

   do iz = 1, nz+1
      zfaces(iz) = (iz - 1)*dz
      H_exact (iz) = sin(zfaces(iz))
      Hp_exact(iz) = cos(zfaces(iz))
   end do

   ! Shu-Osher method:
   call flux_diff(nz, dz, F, Fd, H, Hp)

   open (unit = 1, file = 'Fd_Shu_Osher.dat',  form = 'formatted', status = 'unknown')
   open (unit = 2, file = 'err_Fp_Shu_Osher.dat', form = 'formatted', status = 'unknown')
   open (unit = 3, file = 'err_H.dat',  form = 'formatted', status = 'unknown')
   open (unit = 4, file = 'err_Hp.dat', form = 'formatted', status = 'unknown')
   open (unit = 7, file = 'H.dat',      form = 'formatted', status = 'unknown')
   open (unit = 8, file = 'Hp.dat',     form = 'formatted', status = 'unknown')

   do iz = 1, nz
      write(1, "(3(1x, e12.5))") zcen  (iz), F(iz), Fd(iz)
      write(2, "(2(1x, e12.5))") zcen  (iz), Fd(iz) - exact_deriv(iz)
   end do
   print *, ' Shu-Oshe error at iz = 1 is ',  Fd(1) - exact_deriv(1)
   write (9, "(2(1x, e12.5))") FLOAT(nz), ABS(Fd(1) - exact_deriv(1))

   do iz = 1, nz+1
      write(3, "(2(1x, e12.5))") zfaces(iz), H(iz ) - H_exact (iz)
      write(4, "(2(1x, e12.5))") zfaces(iz), Hp(iz) - Hp_exact(iz)
      write(7, "(3(1x, e12.5))") zfaces(iz), H (iz), H_exact (iz)
      write(8, "(3(1x, e12.5))") zfaces(iz), Hp(iz), Hp_exact(iz)
   end do
   write (11, "(2(1x, e12.5))") FLOAT(nz), ABS(Hp(1) - Hp_exact(1))

   close(1)
   close(2)
   close(3)
   close(4)
   close(7)
   close(8)

   ! Straight Pade differentiation (non_Shu):
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call pade_diff(nz, dz, F, Fd)

   open (unit = 1, file = 'err_Fp_straight_Pade.dat', form = 'formatted', status = 'unknown')

   do iz = 1, nz
      write(1, "(3(1x, e12.5))") zcen(iz), ABS(Fd(iz) - exact_deriv(iz))
   end do
   print *, ' Straight Pade error at iz = 1 is ', Fd(1) - exact_deriv(1)
   write (10, "(2(1x, e12.5))") FLOAT(nz), ABS(Fd(1) - exact_deriv(1))

   close(1)
   print *, ' wrote err_straight_Pade.dat with nz = ', nz

   ! Differentiation of H_exact in main program:
   call pade_diff(nz+1, dz, H_exact, Hp)
   write (12, "(2(1x, e12.5))") FLOAT(nz), ABS(Hp(1) - Hp_exact(1))

   ! Differentiation of sin at the cell centers:

   deallocate (zcen, F, Fd, exact_deriv, cc_nz)
   deallocate (zfaces, H, Hp, H_exact, Hp_exact, cc_nzp1)
end do

close(9)
close(10)
close(11)
close(12)

stop
end

!----------------------------------------------------------------------------------85

subroutine flux_diff(n, dz, F, Fd, H, Hp)

! Given the flux at the midpoint of cells, obtains the a flux difference at the same
! point.  This flux difference is a telescoping sum using the trick of Shu and Osher.
! Oct. 24, 2017.

implicit none
integer, intent(in) :: n ! number of cells
real(8), intent(in) :: dz
real(8), dimension(n), intent(in)  :: F  ! Fluxes at cell midpoints.
real(8), dimension(n), intent(out) :: Fd ! Flux differences at the same points.
real(8), dimension(n+1), intent(out) :: H, Hp ! For testing.

! Locals:
integer :: i

! Compute the quantity "H" in the Shu and Osher trick, which is simply the accumulated
! flux.  Computing it seems to be unavoidable.

H(1) = 0.0d0
do i = 2, n+1
   H(i) = H(i-1) + dz*F(i-1)
end do

call pade_diff(n+1, dz, H, Hp)

do i = 1, n
   Fd(i) = (Hp(i+1) - Hp(i)) / dz
end do

return
end

!----------------------------------------------------------------------------------85

subroutine pade_diff(n, dx, f, fp)

! Standard 4th order Pade differentiation with 4th order numerical boundary scheme.
! Oct. 24, 2017.

implicit none
integer,               intent(in)  :: n
real(8)                            :: dx
real(8), dimension(n), intent(in)  :: f
real(8), dimension(n), intent(out) :: fp

! Coefficients for interior scheme:
real(8), parameter :: alpha = 0.25d0, a_over_2 = 0.75d0

! Coefficients for boundary scheme (4.1.4 in Lele's paper) and using
! Lele's notation:
real(8), parameter :: beta = 3.d0, a1 = -17.d0/6.d0, b1 = 1.5d0, c1 = 1.5d0, &
     d1 = -1.d0/6.d0

! Locals:
integer :: i

! This is the RHS which will eventually become the the derivative.
fp(1) =  a1*f(1) + b1*f(2  ) + c1*f(3  ) + d1*f(4  )
fp(n) = -a1*f(n) - b1*f(n-1) - c1*f(n-2) - d1*f(n-3)

do i = 2, n - 1
   fp(i) = a_over_2 * (f(i+1) - f(i-1))
end do

! Call is:
!tridiag (m, a, b1, c1, bm, am, x)
! Don't confuse with Lele's notation.
! First row is b1, c1 = 1, beta
! Last row is am, bm = beta 1

call tridiag(n, alpha, 1.0d0, beta, 1.0d0, beta, fp)

do i = 1, n
   fp(i) = fp(i)/dx
end do

return
end subroutine pade_diff

!----------------------------------------------------------------------------------85

subroutine tridiag(m, a, b1, c1, bm, am, x)

! Tridiagonal solver.  Courtesy of Alan Wray's routine triddx.

implicit none
integer, intent(in)    :: m ! # of rows
real(8), intent(in)    :: a, b1, c1, bm, am
real(8), intent(inout) ::  x(m)

! Local:
integer :: j
real(8) :: cc(m)

  !   Solves the following special tridiagonal system; r is passed in x
  !
  !    (b1 c1 0 . . . . . . . ) x(1)      r(1)
  !         
  !    (a  1  a  0 . . . . .  ) x(2)      r(2)
  !           
  !    (0  a  1  a  0 . . . . ) x(3)  =   r(3)
  !            
  !     . . . . . . . . . .
  !
  !    (. . . .  0   a   1  a ) x(m-1)    r(m-1)
  !                 
  !    (. . . . . .  0   am bm) x(m)      r(m)

cc(1) = c1/b1
do j = 2, m-1
   cc(j) = a/(1.d0 - a*cc(j-1))
end do

x(1) = x(1)/b1

do j = 2, m-1
   x(j) = (x(j) - a*x(j-1))/(1.d0 - a*cc(j-1))
end do

x(m) = (x(m) - am*x(m-1))/(bm - am*cc(m-1))

Do j = m-1, 1, -1
   x(j) = x(j) - cc(j)*x(j+1)
end Do

end subroutine tridiag

!----------------------------------------------------------------------------------85
