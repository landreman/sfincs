# Makefile for a library that reads VMEC wout files.
# Based on the JMC code by Joachim Geiger at IPP-Greifswald.
# Uses PPPL's ezcdf routines (http://w3.pppl.gov/ntcc/EZcdf).
# Typically this makefile is only called from make in the parent directory,
# in which the variables FC and EXTRA_COMPILE_FLAGS are set.

OBJ_FILES = \
	ezcdf_attrib.o \
	ezcdf.o \
	ezcdf_genget.o \
	ezcdf_genput.o \
	ezcdf_inqvar.o \
	ezcdf_opncls.o \
	handle_err.o \
	constants.o \
	kind_defs.o \
	pisa_io.o \
	pisa_vmec_module.o \
	readVMEC.o

all: readVMEC.a

readVMEC.a: $(OBJ_FILES)
	ar rcs readVMEC.a $(OBJ_FILES)
#	rm *.o

#all: $(OBJ_FILES)

%.o: %.f90
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
pisa_io.o: kind_defs.o
pisa_vmec_module.o: kind_defs.o pisa_io.o ezcdf.o
#wout_data.o: kind_defs.o constants.o pisa_vmec_module.o
readVMEC.o: pisa_vmec_module.o
