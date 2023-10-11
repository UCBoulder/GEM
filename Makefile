SRCS =	gem_com.f90 gem_equil.f90 gem_main.f90 gem_outd.f90 gem_fcnt.f90 gem_fft_wrapper.f90 gem_gkps_adi.f90 gem_beam.f90

OBJS =	gem_com.o gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_fft_wrapper.o gem_gkps_adi.o gem_beam.o

#LIBS = ~/installed/dfftpack_cray/libdfftpack.a -lfftw3_threads -lfftw3  -lfftw3f_threads -lfftw3f
LIBS = /global/homes/u/u10198/installed/dfftpack_cray/libdfftpack.a
PLIB = gem_pputil.o


# -acc=gpu by default
# -ta is deprecated
#OPT = -O0 -r8 -Kieee -llapack -lblas -cpp -Mbounds -g -Minfo=acc -acc -ta=nvidia:cc80
#OPT = -O0 -r8 -Kieee -llapack -lblas -cpp -Mbounds -g -mp #-Minfo=accel
OPENMP_OPT := -DOPENMP
OPENACC_OPT := -DOPENACC

#F90 = ftn $(OPENACC_OPT)#-DOLD_PMOVE #$(OPENACC_OPT)
F90 = ftn $(OPENMP_OPT)#-DOLD_PMOVE #$(OPENMP_OPT)

OPT =  -f free -s real64  -h omp -hvector0 -h acc -hlist=a -e Z
#OPT = -C -g -traceback -Mbounds -r8 -Kieee -llapack -cpp -mp # -lblas -Minfo=accel -mp # -acc
#OPT = -O3 -r8 -Kieee -llapack -cpp -mp # -lblas -Minfo=accel -mp # -acc

#OPT = -O0 -r8 -Kieee -llapack -lblas -cpp -Minfo=accel -mp
#OPT = -O0 -r8 -Kieee -llapack -lblas -cpp -Minfo=accel -acc=multicore
#OPT = -O0 -r8 -Kieee -llapack -lblas -cpp -Minfo=accel -acc=host

LDFLAGS = 

#all : gem

gem_main: gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_pputil.o gem_com.o gem_fft_wrapper.o gem_gkps_adi.o gem_beam.o
	$(F90)  -o gem_main $(OPT) $(OBJS) $(PLIB) $(LIBS) 

gem_pputil.o: gem_pputil.f90
	$(F90) -c $(OPT) gem_pputil.f90

gem_com.o: gem_com.f90 gem_pputil.o
	$(F90) -c $(OPT) gem_com.f90

gem_equil.o: gem_equil.f90 gem_pputil.o gem_com.o
	$(F90) -c $(OPT) gem_equil.f90

gem_gkps_adi.o: gem_gkps_adi.f90 gem_com.f90 gem_equil.f90 gem_pputil.f90
	$(F90) -c $(OPT) gem_gkps_adi.f90

gem_main.o: gem_main.f90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o gem_gkps_adi.o gem_beam.o
	$(F90) -c $(OPT) gem_main.f90

gem_outd.o: gem_outd.f90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o
	$(F90) -c $(OPT) gem_outd.f90

gem_fcnt.o: gem_fcnt.f90
	$(F90) -c $(OPT) gem_fcnt.f90

gem_fft_wrapper.o: gem_fft_wrapper.f90
	$(F90) -c $(OPT) gem_fft_wrapper.f90

gem_beam.o: gem_beam.f90 gem_pputil.o gem_com.o gem_equil.o
	$(F90) -c $(OPT) gem_beam.f90

clean:
	rm -f *.o *.lst *.mod gem_main
