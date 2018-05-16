# Common commands
# -------------------------------
# To build SFINCS:
#   make
# or, equivalently,
#   make all
# To delete .mod and .o files:
#   make clean

# Do you want to link to NetCDF?  Doing so is necessary to read VMEC wout files in NetCDF
# format, but not necessary to read VMEC wout files in ASCII format, or for other geometry input.
# If you do want to link NetCDF, set the USE_NETCDF variable here to any value; the actual value
# does not matter.  If you do not want to link NetCDF, do not assign any value to USE_NETCDF here.
USE_NETCDF ?= X

##############################################################################3
# End of parameters
##############################################################################3

# First, make sure the environment variable SFINCS_SYSTEM is set appropriately:

ifndef SFINCS_SYSTEM
$(error You must set the environment variable SFINCS_SYSTEM to select which makefile to use. For example, to use makefiles/makefile.edison, set SFINCS_SYSTEM=edison. To set this environment variable, use the export or setenv command, depending on your shell)
endif

SYSTEM_MAKEFILE = ../makefiles/makefile.$(SFINCS_SYSTEM)
SYSTEM_SFINCS_MAIN = ..

ifeq ($(wildcard $(SYSTEM_MAKEFILE)),)
$(error Based on the environment variable SFINCS_SYSTEM, make is trying to include the file $(SYSTEM_MAKEFILE), but this file does not exist)
endif

.PHONY: all clean test retest

all: sfincs_testAdjointDiagnostics

include $(SYSTEM_MAKEFILE) 

ifdef USE_NETCDF
	EXTRA_COMPILE_FLAGS += -DUSE_NETCDF
endif

include makefile.depend

LIBSTELL_DIR = ../mini_libstell
LIBSTELL_FOR_SFINCS = ../mini_libstell/mini_libstell.a
SFINCS_LIB = ../libsfincs.a
# Export variables that will be used in the mini_libstell subdirectory or for "make test":
export

%.o: %.F90
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -I$(SYSTEM_SFINCS_MAIN) -c $<

%.o: %.f
	$(FC) $(EXTRA_COMPILE_FLAGS) -c $<

sfincs_testAdjointDiagnostics: sfincs_testAdjointDiagnostics.o
	$(FLINKER) -o sfincs_testAdjointDiagnostics sfincs_testAdjointDiagnostics.o $(LIBSTELL_FOR_SFINCS) $(SFINCS_LIB) $(EXTRA_LINK_FLAGS) $(PETSC_LIB)

# Since the sometimes-included PETSc makefile uses a double-colon to define "clean", we must do the same here:
clean::
	rm -f *.o *.mod *.MOD *.a *.dat *~ sfincs_testAdjointDiagnostics
