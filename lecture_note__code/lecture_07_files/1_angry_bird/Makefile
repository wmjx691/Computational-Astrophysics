F90=gfortran
F90FLAGS= -O2 -fdefault-real-8 -fdefault-integer-8

SOURCES=angry.f90 constants.f90 physics.f90 solver.f90

CLEANSTUFF=rm -rf *.o *.mod *~ angry

OBJECTS=$(SOURCES:.f90=.o)

all: mod angry

angry: $(OBJECTS)
	$(F90) $(F90FLAGS) $(OBJECTS) -o angry

mod: constants.f90 solver.f90 physics.f90
	$(F90) $(F90FLAGS) -c $^

$(OBJECTS): %.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	$(CLEANSTUFF)
