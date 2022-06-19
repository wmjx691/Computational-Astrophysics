#
#
#

F90=gfortran
F90FLAGS= -O2 -fdefault-real-8 -fdefault-integer-8

#
#
#
SOURCES= Simulation_data.f90 IO.f90 initial.f90 boundary.f90 finalize.f90 evolution.f90 advection.f90 

CLEANSTUFF=rm -rf *.o *.mod *~ *.d advection

OBJECTS=$(SOURCES:.f90=.o)

all: advection

advection: $(OBJECTS)
	$(F90) $(F90FLAGS) $(OBJECTS) -o advection

$(OBJECTS): %.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	$(CLEANSTUFF)
