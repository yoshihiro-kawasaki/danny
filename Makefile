FC = gfortran

# FLAG = -c -O3 -fallow-argument-mismatch

FLAG = -c -O3 

OBJDIR = ./obj

OBJ = src/physical_constants.o \
	  src/danny_variables.o \
	  src/danny_read.o \
	  src/danny_dust.o \
	  src/danny.o \
	  src/main.o \
	  src/odepack/opkdmain.o \
	  src/odepack/opkda1.o \
	  src/odepack/opkda2.o 

run: $(OBJ)
	$(FC) -o run $(OBJ)

src/physical_constants.o: src/physical_constants.f90
	(cd src; $(FC) $(FLAG) physical_constants.f90)

src/danny_variables.o: src/danny_variables.f90
	(cd src; $(FC) $(FLAG) danny_variables.f90)

src/danny_read.o: src/danny_read.f90 
	(cd src; $(FC) $(FLAG) danny_read.f90)

src/danny.o: src/danny.f90 
	(cd src; $(FC) $(FLAG) danny.f90)

src/danny_dust.o: src/danny_dust.f90 
	(cd src; $(FC) $(FLAG) danny_dust.f90)
	
src/main.o: src/main.f90 
	(cd src; $(FC) $(FLAG) main.f90)

src/odepack/opkdmain.o: src/odepack/opkdmain.f
	(cd src/odepack; $(FC) $(FLAG) opkdmain.f)

src/odepack/opkda1.o: src/odepack/opkda1.f
	(cd src/odepack; $(FC) $(FLAG) opkda1.f)

src/odepack/opkda2.o: src/odepack/opkda2.f 
	(cd src/odepack; $(FC) $(FLAG) opkda2.f)

# src/odepack/opkdmain.o: src/odepack/opkdmain.f90 
# 	(cd src/odepack; $(FC) $(FLAG) opkdmain.f90)

# src/odepack/opkda1.o: src/odepack/opkda1.f90 
# 	(cd src/odepack; $(FC) $(FLAG) opkda1.f90)

# src/odepack/opkda2.o: src/odepack/opkda2.f90 
# 	(cd src/odepack; $(FC) $(FLAG) opkda2.f90)

clean:
	rm run $(OBJ)