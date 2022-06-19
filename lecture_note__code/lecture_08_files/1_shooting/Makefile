#
#

F90=gfortran
F90FLAGS= -O2 -fdefault-real-8 -fdefault-integer-8

######

SOURCES=solver.f90 example1.f90 example2.f90
 
CLEANSTUFF=rm -rf *.o *.mod example1 example2 *~

all: solver example1 example2 

solver: Solver.f90
	$(F90) $(F90FLAGS) -c $^

example2: Solver.f90 example2.f90
	$(F90) $(F90FLAGS) Solver.f90 example2.f90 -o example2

example1: Solver.f90 example1.f90
	$(F90) $(F90FLAGS) Solver.f90 example1.f90 -o example1

$(OBJECTS): %.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	$(CLEANSTUFF)
