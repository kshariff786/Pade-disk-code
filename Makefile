# You can type the following to the terminal:
# make pade        ! To make the serial   executable called pade
# make pade_mpi    ! To make the parallel executable called pade_mpi.  You can have objects and executable of both sitting in the directory.
# make clean       ! Deletes the executables (pade and/or pade_mpi) and all .o files.
# make pade mode=debug ! To make the serial executable with debug output and array bounds checking turned on.  This will consume extra CPU time.
# make pade_mpi mode=debug ! To make the parallel executable with debug output and runtime bounds
# checking

# make pade transpose_timing=yes or make pade_mpi transpose_timing=yes:
# This will perform cpu timing for transposes.  The result will be output to
# stdout at the end of every step and at the end of the run.

################################################################################
#                    Code options
phi_averaged_statistics = yes

# I need the math77 library for Bessel function evaluations for the Seligman vortex:
# This is a flag.  Set it equal to yes or no.  If you set it to "no" you will
# need to remove app_vortex_pair.o from being linked.
math77 = yes

################################################################################

# Set this to the shell you use.  I need this because I have gfortran aliased in my .tcshrc
SHELL = /bin/tcsh

# In the next two lines, specify how to invoke the serial and parallel fortran
# compilers on your machine.  If you only want to make the serial code (pade)
# then don't worry about the mpi_fortran variable.

ifeq ($(OSTYPE),darwin)
#  For macs:
   serial_fortran = gfortran
   mpi_fortran = mpif90
else ifeq ($(HOST),linux261.nas.nasa.gov)
#  For my linux box:
   serial_fortran = gfortran
   mpi_fortran = mpif90
else ifeq ($(MYHOST),pfe)
# Note: In my .cshrc I have setenv MYHOST pfe
#  For the Pleiades supercomputer at NASA Ames:
   serial_fortran = ifort
   mpi_fortran = mpif90
endif

# In the next line, select the fft routine to be used by the corrected FARGO trick.
# The two choices are fftw or rogallo.  For fftw you need to specify the include paths
# below.
#fft = rogallo
fft = fftw

# Set this for your situation (needed only if you chose fftw above).
ifeq ($(fft),fftw)
   ifeq ($(OSTYPE),darwin)
      fftw_include_path = $(HOME)/Applications/fftw-3.3.7/include
      fftw-library-path = $(HOME)/Applications/fftw-3.3.7/lib
   else ifeq ($(HOST),linux261.nas.nasa.gov)
      fftw_include_path = $(HOME)/WORK/fftw-3.3.8/include
      fftw-library-path = $(HOME)/WORK/fftw-3.3.8/lib
   else
#     Pleiades:
      fftw_include_path = $(HOME)/Applications/fftw-3.3.7/include
      fftw-library-path = $(HOME)/Applications/fftw-3.3.7/lib
   endif
endif

# If you set math77 to yes then you need to supply the location of the library.
# Note: Pleiades has the math77 library on the system so you don't need to do anything.
ifeq ($(math77),yes)
   ifeq ($(OSTYPE),darwin)
      math77-library-path = $(HOME)/Applications/MATH77
   else ifeq ($(HOST),linux261.nas.nasa.gov)
      math77-library-path = $(HOME)/WORK/MATH77
   endif
   library_links += -L$(math77-library-path) -lmath77
endif

# The targets for this makefile are:
#    clean : Delete all objects, executables, and compiled modules.

#    pade : Serial pade code
#    pade_mpi : Parallel pade code
#
#    To do: Remove these
#    flux_diff_test : Serial code to test the flux differentiation method.
#    Shu_Osher_test : Serial code to test the Shu-Osher method for getting
#                     a telescoping sum derivative.
#    phi_advection_test : An incomplete code

# To compile in debug mode type for instance:
# make pade mode=debug
# Otherwise production mode will be assumed.  In debug mode you will get
# output about what is going on and compiler run-time checks are turned-on.
# This will slow things down.
#
# You don't have to read any further if you don't want to.

# The extension .mpio is used for object files pertaining to the mpi code.
# This is to allow both the serial and mpi codes to be built in the same
# directory.
.SUFFIXES : .f .f90 .o .mpio

ifeq ($(MAKECMDGOALS),pade_mpi)
   fortran = $(mpi_fortran)
else
   fortran = $(serial_fortran)
endif

# List of preprocessor options we sometimes use:
enable_preprocessing = -cpp
debug_print = -Ddebug_print

# Run time check for debug mode.
ifeq ($(serial_fortran),gfortran)
   debug_mode_run_time_checks = -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow
else ifeq ($(serial_fortran),ifort)
   debug_mode_run_time_checks = -check -g -traceback
endif

# List of possible compiler options I like to use for gfortran:
init_nan = -finit-real=nan
warn_un_init = -Wuninitialized
# Warn about array temporaries created by the compiler
warn_array_temps = -Warray-temporaries
# This has the effect of putting implicit none at the start of every subroutine
# in case I forget.
implicit_none = -fimplicit-none
# In case needed for linking with c objects.
no_underscoreing = -fno-underscoring -fno-second-underscoring
max_errors = -fmax-errors=10

# Directory in which .mod files will be kept.  I simply do not like to see all those
# .mod files mixed in with the .f90 files.
ifeq ($(serial_fortran),gfortran)
   compiled_modules_dir = -JCOMPILED_MODULES
else ifeq ($(fortran),ifort)
   compiled_modules_dir = -module COMPILED_MODULES
endif

optimization_level = 3

# For sanity checks, the code checks for uninitialized variables by comparing against NaN, so
# DO NOT remove the init_nan or init=snan option below:

# Fortran options I always use:
ifeq ($(serial_fortran),gfortran)
   fortran_options = $(init_nan) $(warn_un_init) $(implicit_none) $(compiled_modules_dir) $(enable_preprocessing) $(warn_array_temps) $(max_errors)
else ifeq ($(serial_fortran),ifort)
   fortran_options = -init=snan $(enable_preprocessing) $(compiled_modules_dir)
endif

# Add some options for debug mode:
ifeq ($(mode), debug)
   fortran_options += $(debug_print) $(debug_mode_run_time_checks)
endif

# Add pre-processor variable if needed:
ifeq ($(transpose_timing), yes)
   fortran_options += -Dtranspose_timing
endif


# Add options for Pleiades to optimize for the processors used:
#ifeq ($(MYHOST),pfe)
#   fortran_options += -axCORE-AVX2 -xSSE4.2 -ip
#endif

! For FFTW to implement fargo shifts:
ifeq ($(fft), fftw)
   fortran_options += -I$(fftw_include_path)
   library_links += -L$(fftw-library-path) -lfftw3 -lm
#  This adds an ifdef preprocessor option called fftw for the code to use:
   fortran_options += -Dfftw
endif

ifeq ($(fft), fftw)
ifeq ($(MAKECMDGOALS),pade_mpi)
   library_links += -L$(fftw-library-path) -lfftw3_mpi
endif
endif

ifeq ($(MAKECMDGOALS), pade_mpi)
   fortran_options += -Dmpi_code
endif

# If we are using the ifort compiler then the code needs to issue a
# use IFPORT
# the ifort portability library in order than the gfortran random number
# generator can be used.  This pre-processor variable allows us to do that.
ifeq ($(serial_fortran), ifort)
   fortran_options += -Difort
endif

exec1 = pade
exec2 = pade_mpi

# These were used for micro testing purposes while the code was being developed.
exec3 = flux_diff_test
exec4 = Shu_Osher_test
exec5 = phi_advection_test

# The order of the .o files is important to ensure that the modules required by succeeding
# objects have been compiled.

# Objects for serial pade code.  The main program is in pade.f90
objects_for_pade =	modules.o \
			pade.o \
			basic_state.o \
			smooth.o \
			rogallo_fft.o \
			fargo_trick_and_plotting_shift.o \
			boundary_conditions.o \
			viscous_terms.o \
			artificial_pressure.o \
			activate_routines.o \
			set_up_routines.o \
			sponge.o \
			gravity.o \
			pade_diff.o \
			grid.o \
			domain_mesh_and_partition.o \
			add_derivative_routines.o \
			rotating_frame.o \
			rhs.o \
			control_and_check_routines.o \
			time_stepping_routines.o \
			tridiagonal_solvers.o \
			transpose_routines.o \
			plotting_output.o \
			restart_and_save_file_operations.o \
			vorticity_and_dilatation.o \
			baroclinic_term.o \
			Cassen_Moosman_infall.o \
			zeroin.o \
			conservation_diagnostics.o

# The appropriate subroutine is invoked in the main program in
# source file pade.f90 according to the i_run_type variable
# specified in the input file.
application_objects =	app_user.o \
			app_euler1d_tests.o \
			app_hydrostatic_test.o \
			app_homentropic_solid_body_rotation_test.o \
			app_vsi_3D.o \
			app_vertically_integrated_disk.o \
			app_vortex_pair.o \
			app_taylor_couette.o \
			app_zvi.o

objects_for_pade += $(application_objects)

ifeq ($(phi_averaged_statistics), yes)
   objects_for_pade += phi_avg_statistics.o
endif

# The objects for the mpi code are the same except that they have the extension .mpio
objects_for_pade_mpi = $(objects_for_pade:.o=.mpio)

# For micro testing purposes as the code was being developed.
objects3 = flux_diff_test.o flux_diff.o tridiagonal_solvers.o
objects4 = Shu_Osher_test.o
objects5 = phi_advection_test.o flux_diff.o tridiagonal_solvers.o

#

$(exec1):	$(objects_for_pade)
		$(fortran) $(compiled_modules_dir) -o $(exec1) $(objects_for_pade) $(library_links)

$(exec2):	$(objects_for_pade_mpi)
		$(mpi_fortran) $(compiled_modules_dir) -o $(exec2) $(objects_for_pade_mpi) $(library_links)

$(exec3):	$(objects3)
		$(fortran) $(compiled_modules_dir) -o $(exec3) $(objects3)

$(exec4):	$(objects4)
		$(fortran) $(compiled_modules_dir) -o $(exec4) $(objects4)

$(exec5):	$(objects5)
		$(fortran) $(compiled_modules_dir) -o $(exec5) $(objects5)

$(exec6):	$(objects6)
		$(fortran) $(compiled_modules_dir) -o $(exec6) $(objects6)

# The dashes below make the command continue even after an error, e.g., if the file does not exist.
clean:
	-\rm COMPILED_MODULES/*.mod
	-\rm *.mod
	-\rm *.o
	-\rm *.mpio
	-\rm $(exec1) $(exec2) $(exec3) $(exec4) $(exec5) $(exec6)

.f90.o :
	echo $(fortran_options)
	$(serial_fortran) -O$(optimization_level) $(fortran_options) -c $<

.f90.mpio :
	$(mpi_fortran) -O$(optimization_level) $(fortran_options) -c -o$(<:.f90=.mpio) $<

.f.mpio :
	$(mpi_fortran) -O$(optimization_level) $(fortran_options) -c -o$(<:.f=.mpio) $<

.f.o :
	$(serial_fortran) -O$(optimization_level) $(fortran_options) -c $<








