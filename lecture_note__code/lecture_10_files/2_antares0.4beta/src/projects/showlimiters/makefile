default:
	ifort -c *.f90
	ifort -c *.f
	ifort *.o -o showlimiters.exe
	mv showlimiters.exe ../bin
clean:	
	rm *.o
	rm *.exe
	rm *.tab
