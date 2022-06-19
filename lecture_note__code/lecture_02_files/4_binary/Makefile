F90=gfortran
F90FLAGS= -O2 -fdefault-real-8 -fdefault-integer-8

SOURCES=binary.f90 constants.f90 physics.f90 Simulation_data.f90 output.f90

CLEANSTUFF=rm -rf *.o *.mod *.dat *~ binary

OBJECTS=$(SOURCES:.f90=.o)

all: mod binary

binary: $(OBJECTS)
	$(F90) $(F90FLAGS) $(OBJECTS) -o binary

mod: Simulation_data.f90 constants.f90 physics.f90 output.f90
	$(F90) $(F90FLAGS) -c $^

$(OBJECTS): %.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	$(CLEANSTUFF)
