#compiler options for debugging
#c=gfortran -fimplicit-none -fbounds-check -fbacktrace -g -O0 -Wunused -finit-real=nan
c=gfortran -fbounds-check -fbacktrace -g -O0 -Wunused -finit-real=nan
#compiler options for fast execution
#c=gfortran -03

all_objs:= types.o functions.o tools.o parameters.o fdm.o

all: $(all_objs) main.o
	$c -g -o HeatFF $(all_objs) main.o


#all: types.o functions.o tools.o main.o # dependencies
#	$c -g -o gauss_roundoff types.o functions.o tools.o main.o # rules how the file is created

types.o: types.f90
	$c -c types.f90

functions.o: types.o functions.f90
	$c -c functions.f90

tools.o: tools.f90
	$c -c tools.f90

parameters.o: types.o parameters.f90
	$c -c parameters.f90

fdm.o: types.o tools.o fdm.f90
	$c -c fdm.f90

main.o: types.o functions.o tools.o main.f90
	$c -c main.f90

clean: 	
	rm -rf *.mod *.o
