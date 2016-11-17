# Makefile for the quadpack library.
# Typically this makefile is only called from make in the parent directory,
# in which the variables FC and EXTRA_COMPILE_FLAGS are set.

OBJ_FILES = d1mach.o  dqagi.o  dqk15i.o  dqk41.o  dqpsrt.o   qage.o   qk15.o   qk31.o  qk61.o \
	dqage.o   dqelg.o  dqk21.o   dqk51.o qagie.o  qk15i.o  qk41.o  qpsrt.o \
	dqagie.o  dqk15.o  dqk31.o   dqk61.o  qelg.o   qk21.o   qk51.o  r1mach.o


all: $(OBJ_FILES)
	ar rcs quadpack.a $(OBJ_FILES)
	rm *.o

%.o: %.f
	$(FC) $(EXTRA_COMPILE_FLAGS) -c $<

clean::
	rm -f *.o *.a