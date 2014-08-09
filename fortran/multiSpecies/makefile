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

# First, make sure the environment variable SFINCS_SYSTEM is set appropriately:

ifndef SFINCS_SYSTEM
$(error You must set the environment variable SFINCS_SYSTEM to select which makefile to use. For example, to use makefiles/makefile.edison, set SFINCS_SYSTEM=edison. To set this environment variable, use the export or setenv command, depending on your shell)
endif

SYSTEM_MAKEFILE = makefiles/makefile.$(SFINCS_SYSTEM)
ifeq ($(wildcard $(SYSTEM_MAKEFILE)),)
$(error Based on the environment variable SFINCS_SYSTEM, make is trying to include the file $(SYSTEM_MAKEFILE), but this file does not exist)
endif

TARGET = sfincs

OBJ_FILES = sfincs.o xGrid.o writeHDF5Output.o globalVariables.o scan.o readInput.o \
	polynomialDiffMatrices.o uniformDiffMatrices.o \
	uniformInterpolationMatrix.o printToStdout.o  \
	indices.o sparsify.o geometry.o \
	polynomialInterpolationMatrix.o solveDKE.o quadpack/quadpack.a 

.PHONY: all clean test retest

all: $(TARGET)

include $(SYSTEM_MAKEFILE)

# Export variables that will be used in the quadpack subdirectory or for "make test":
#export FC
#export EXTRA_COMPILE_FLAGS
#export SFINCS_RUN_LARGE_EXAMPLES
export

%.o: %.F90
	$(FC) $(EXTRA_COMPILE_FLAGS) -c $<

quadpack/quadpack.a:
	$(MAKE) -C quadpack

$(TARGET): $(OBJ_FILES)
	$(FLINKER) -o $(TARGET) $(OBJ_FILES) $(EXTRA_LINK_FLAGS) $(PETSC_LIB)

# Since the sometimes-included PETSc makefile uses a double-colon to define clean, we must do the same here:
clean::
	rm -f *.o *.mod *.a *.dat *~ $(TARGET) sfincsMatrices.m sfincsOutput.h5
	cd quadpack; rm -f *.o *.a

test_make:
	@echo SFINCS_SYSTEM is $(SFINCS_SYSTEM)
	@echo FC is $(FC)
	@echo FLINKER is $(FLINKER)
	@echo EXTRA_COMPILE_FLAGS is $(EXTRA_COMPILE_FLAGS)
	@echo EXTRA_LINK_FLAGS is $(EXTRA_LINK_FLAGS)
	@echo PETSC_LIB is $(PETSC_LIB)
	@echo TARGET is $(TARGET)
	@echo OBJ_FILES is $(OBJ_FILES)

test: $(TARGET)
	@echo "Beginning functional tests." && cd examples && export SFINCS_RETEST=no && ./runExamples.py

retest: $(TARGET)
	@echo "Testing sfincsOutput.h5 for examples without re-running then." && cd examples && export SFINCS_RETEST=yes && ./runExamples.py

# Dependencies:
writeHDF5Output.o: globalVariables.o scan.o
printToStdout.o: globalVariables.o
readInput.o: globalVariables.o
geometry.o: globalVariables.o
solveDKE.o: globalVariables.o xGrid.o polynomialDiffMatrices.o printToStdout.o sparsify.o geometry.o indices.o
indices.o: globalVariables.o
scan.o: globalVariables.o
sfincs.o: globalVariables.o writeHDF5Output.o printToStdout.o solveDKE.o readInput.o scan.o

