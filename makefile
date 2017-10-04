EXEC = power.e

SYSTYPE="emiliano"

ifeq ($(SYSTYPE),"emiliano")
FC          =  gfortran
FLAGS       = -O3 -Wall 
FFTW_LIBR   = -L/usr/local/lib -lfftw3 -lm
FFTW_INCL   = -I/usr/local/include #-I/Users/Emiliano/Cosmo/Tools/FFTW/fftw-3.3.4/api/
endif

ifeq ($(SYSTYPE),"martin")
FC          =  ifort  
FLAGS       = -m64 
FFTW_LIBR   = -L/Users/martincrocce/local/lib -lfftw3 -lm
FFTW_INCL   = -I/Users/martincrocce/local/include/        # path to fftw3.f
endif

MODULES =  parbox.o grid.o interpolation.o input_catalog.o #fftw3.o
MAIN = PowerI4.o

%.o: %.f90
	$(FC) $(FLAGS)  $(FFTW_INCL) -c $<

power: $(MODULES) $(MAIN)
	$(FC) $(FLAGS) $(FFTW_INCL) -o $(EXEC) $(MODULES) $(MAIN) $(FFTW_LIBR)

clean:
	\rm -f $(EXEC) *.o *.mod

# makefile by Geppetto  
