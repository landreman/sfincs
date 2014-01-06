# Makefile for the quadpack library

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

OBJS = d1mach.o  dqagi.o  dqk15i.o  dqk41.o  dqpsrt.o   qage.o   qk15.o   qk31.o  qk61.o \
	dqage.o   dqelg.o  dqk21.o   dqk51.o qagie.o  qk15i.o  qk41.o  qpsrt.o \
	dqagie.o  dqk15.o  dqk31.o   dqk61.o  qelg.o   qk21.o   qk51.o  r1mach.o


%.o: %.f
	$(FLINKER) -c $<

all: $(OBJS)
	ar rcs quadpack.a $(OBJS)
	rm *.o

cleaner:
	rm -f *.o *.a