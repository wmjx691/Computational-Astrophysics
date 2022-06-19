F90=gfortran
F90FLAGS= -O2 -fdefault-real-8 -fdefault-integer-8

SOURCES=matrix.f90 linalg.f90

CLEANSTUFF=rm -rf *.o *.mod *.dat *~ matrix

OBJECTS=$(SOURCES:.f90=.o)

all: mod main

main: $(OBJECTS)
	$(F90) $(F90FLAGS) $(OBJECTS) -o matrix

mod: linalg.f90 
	$(F90) $(F90FLAGS) -c $^

$(OBJECTS): %.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	$(CLEANSTUFF)
