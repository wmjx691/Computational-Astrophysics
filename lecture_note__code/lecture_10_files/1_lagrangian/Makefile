F90=gfortran
F90FLAGS= -O2 -fdefault-real-8 -fdefault-integer-8

SOURCES= Simulation_data.f90 IO.f90 lagrangian.f90

CLEANSTUFF=rm -rf *.o *.mod *~ *.d lagrangian

OBJECTS=$(SOURCES:.f90=.o)

all: lagrangian

lagrangian: $(OBJECTS)
	$(F90) $(F90FLAGS) $(OBJECTS) -o lagrangian

$(OBJECTS): %.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	$(CLEANSTUFF)
