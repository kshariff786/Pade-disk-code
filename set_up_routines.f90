!----------------------------------------------------------------------------------85

subroutine set_up_thermal_parameters(gamma_arg, isothermal_arg)

use grid
use thermal_parameters
use logical_units
use partition_data
use control_parameters
implicit none
real(8) :: gamma_arg
logical :: isothermal_arg

gamma           = gamma_arg
isothermal      = isothermal_arg

gm1 = gamma - 1

! Note: This is for the whole field.  I am having each processor have
! all of r.  I am not sure why I did this.  I'll try to fix it later.
! May be you did it so it could be used in all spaces.
! Note: The user application needs to set this.
allocate (ci_squared_initial(nr, nz))
n_words_allocated = n_words_allocated + nr*nz

! This is for boundary conditions.
if (isothermal) then
   allocate(d_ci_dr_inner(nz), d_ci_dr_outer(nz))
   n_words_allocated = n_words_allocated + 2*nz
end if

#ifdef debug_print
   if (my_node .eq. 0) then
      print *, ' my_node = 0: Finished allocating ci_squared_initial'
      print *, ' # Bytes = ', n_words_allocated/1e6*8, ' Mb'
   end if
#endif

if (my_node .eq. 0) then
   open(unit = lun_info, file = 'info.txt', form = 'formatted', &
        status = 'unknown', access = 'append')
   write(lun_info, *) ' In subroutine set_up_thermal_parameters:'
   write(lun_info, *) ' isothermal = ', isothermal
   close(lun_info)
end if   

end subroutine set_up_thermal_parameters

!----------------------------------------------------------------------------------85
