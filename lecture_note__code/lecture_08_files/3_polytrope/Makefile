#
#

F90=gfortran
F90FLAGS= -O2 -fdefault-real-8 -fdefault-integer-8

######

SOURCES=solver.f90 polytrope.f90 polytrope2.f90
 
CLEANSTUFF=rm -rf *.o *.mod polytrope polytrope2 *~

all: solver polytrope polytrope2

solver: Solver.f90
	$(F90) $(F90FLAGS) -c $^

polytrope2: Solver.f90 polytrope2.f90
	$(F90) $(F90FLAGS) Solver.f90 polytrope2.f90 -o polytrope2

polytrope: Solver.f90 polytrope.f90
	$(F90) $(F90FLAGS) Solver.f90 polytrope.f90 -o polytrope

$(OBJECTS): %.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@




clean:
	$(CLEANSTUFF)
