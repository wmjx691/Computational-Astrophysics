F90=gfortran
F90FLAGS= -O2 -fdefault-real-8 -fdefault-integer-8

SOURCES=nbody.f90 constants.f90 physics.f90 Simulation_data.f90 io.f90

CLEANSTUFF=rm -rf *.o *.mod *.dat *~ nbody

OBJECTS=$(SOURCES:.f90=.o)

all: mod nbody

nbody: $(OBJECTS)
	$(F90) $(F90FLAGS) $(OBJECTS) -o nbody

mod: Simulation_data.f90 constants.f90 io.f90 physics.f90
	$(F90) $(F90FLAGS) -c $^

$(OBJECTS): %.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	$(CLEANSTUFF)
