# Makefile for a library that reads VMEC wout files.
# Based on the JMC code by Joachim Geiger at IPP-Greifswald.
# Uses PPPL's ezcdf routines (http://w3.pppl.gov/ntcc/EZcdf).
# Typically this makefile is only called from make in the parent directory,
# in which the variables FC and EXTRA_COMPILE_FLAGS are set.

OBJ_FILES = \
	constants.o \
	kind_defs.o \
	NEMEC_compute_missing_fields.o \
	pisa_io.o \
	pisa_vmec_module.o \
	read_NEMEC_file.o \
	readVMEC.o

DEPENDENCIES = kind_defs.o pisa_io.o

ifdef USE_NETCDF
	OBJ_FILES += \
	ezcdf.o \
	ezcdf_attrib.o \
	ezcdf_genget.o \
	ezcdf_genput.o \
	ezcdf_inqvar.o \
	ezcdf_opncls.o \
	handle_err.o

	DEPENDENCIES += ezcdf.o
endif

all: readVMEC.a

readVMEC.a: $(OBJ_FILES)
	ar rcs readVMEC.a $(OBJ_FILES)
#	rm *.o

#all: $(OBJ_FILES)

%.o: %.f90
	$(FC) $(EXTRA_COMPILE_FLAGS) -c $<

%.o: %.F90
	$(FC) $(EXTRA_COMPILE_FLAGS) -c $<

clean::
	rm -f *.o *.mod *.a

# Dependencies:
ezcdf_attrib.o: ezcdf_inqvar.o
ezcdf.o: ezcdf_genput.o ezcdf_genget.o ezcdf_attrib.o ezcdf_opncls.o
ezcdf_genget.o: ezcdf_opncls.o ezcdf_inqvar.o
ezcdf_genput.o: ezcdf_opncls.o ezcdf_inqvar.o
ezcdf_inqvar.o:
ezcdf_opncls.o:
handle_err.o:
constants.o: kind_defs.o
kind_defs.o:
NEMEC_compute_missing_fields.o: constants.o kind_defs.o pisa_vmec_module.o
pisa_io.o: kind_defs.o
pisa_vmec_module.o: $(DEPENDENCIES)
read_NEMEC_file.o: constants.o kind_defs.o pisa_vmec_module.o
readVMEC.o: pisa_vmec_module.o
