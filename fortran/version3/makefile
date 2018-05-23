# Common commands
# -------------------------------
# To build SFINCS:
#   make
# or, equivalently,
#   make all
#
# To run tests:
#   make test
#
# To run tests on sfincsOutput.h5 for all examples that would be run for this system with "make test", without actually re-running these examples:
#   make retest
#
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

SYSTEM_MAKEFILE = makefiles/makefile.$(SFINCS_SYSTEM)
ifeq ($(wildcard $(SYSTEM_MAKEFILE)),)
$(error Based on the environment variable SFINCS_SYSTEM, make is trying to include the file $(SYSTEM_MAKEFILE), but this file does not exist)
endif

.PHONY: all clean test retest

all: sfincs 

include $(SYSTEM_MAKEFILE)

ifdef USE_NETCDF
	EXTRA_COMPILE_FLAGS += -DUSE_NETCDF
endif

# DEBUG
EXTRA_COMPILE_FLAGS += -g

include makefile.depend

# Export variables that will be used in the mini_libstell subdirectory or for "make test":
export

%.o: %.F90
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

%.o: %.f
	$(FC) $(EXTRA_COMPILE_FLAGS) -c $<

mini_libstell/mini_libstell.a:
	$(MAKE) -C mini_libstell

# Automatic variables like $@ are documented here: http://www.gnu.org/software/make/manual/make.html#Automatic-Variables
libsfincs.a: $(OBJ_FILES)
	ar rcs $@ $(OBJ_FILES)
#	ar rcs $@ $(OBJ_FILES) mini_libstell/*.o
#	ranlib $@

sfincs: sfincs.o libsfincs.a
	$(FLINKER) -o sfincs sfincs.o libsfincs.a $(LIBSTELL_FOR_SFINCS) $(EXTRA_LINK_FLAGS) $(PETSC_LIB)

# Since the sometimes-included PETSc makefile uses a double-colon to define "clean", we must do the same here:
clean::
	rm -f *.o *.mod *.MOD *.a *.dat *~ sfincs sfincsMatrices.m sfincsOutput.h5
	cd mini_libstell; rm -f *.o *.mod *.MOD *.a
	cd adjointTesting; rm -f *.o *.mod *.a

test_make:
	@echo SFINCS_SYSTEM is $(SFINCS_SYSTEM)
	@echo FC is $(FC)
	@echo CC is $(CC)
	@echo CPP is $(CPP)
	@echo F77 is $(F77)
	@echo CFLAGS is $(CFLAGS)
	@echo CPPFLAGS is $(CPPFLAGS)
	@echo FFLAGS is $(FFLAGS)
	@echo FCFLAGS is $(FCFLAGS)
	@echo LDFLAGS is $(LDFLAGS)
	@echo FLINKER is $(FLINKER)
	@echo EXTRA_COMPILE_FLAGS is $(EXTRA_COMPILE_FLAGS)
	@echo EXTRA_LINK_FLAGS is $(EXTRA_LINK_FLAGS)
	@echo PETSC_LIB is $(PETSC_LIB)
	@echo OBJ_FILES is $(OBJ_FILES)
	@echo LIBSTELL_DIR is $(LIBSTELL_DIR)
	@echo LIBSTELL_FOR_SFINCS is $(LIBSTELL_FOR_SFINCS)

test: sfincs
	@echo "Beginning functional tests." && cd examples && export SFINCS_RETEST=no && ./runExamples.py

retest: sfincs
	@echo "Testing sfincsOutput.h5 for examples without re-running then." && cd examples && export SFINCS_RETEST=yes && ./runExamples.py
