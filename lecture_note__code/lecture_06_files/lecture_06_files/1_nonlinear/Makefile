F90=gfortran
F90FLAGS= -O2 -fdefault-real-8 -fdefault-integer-8

SOURCES=solver.f90 nonlinear.f90

CLEANSTUFF=rm -rf *.o *.mod *.dat *~ nonlinear

OBJECTS=$(SOURCES:.f90=.o)

all: mod main

main: $(OBJECTS)
	$(F90) $(F90FLAGS) $(OBJECTS) -o nonlinear

mod: solver.f90 
	$(F90) $(F90FLAGS) -c $^

$(OBJECTS): %.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	$(CLEANSTUFF)
