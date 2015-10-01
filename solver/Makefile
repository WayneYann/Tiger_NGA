include ../Makefile.in

F90FILES = pressure.f90 hypre.f90 hypre_amg.f90 velocity.f90 \
	   implicit.f90 bodyforce.f90 scalar.f90 scalar_quick.f90 \
	   scalar_weno3.f90 scalar_weno5.f90 fourier.f90 \
	   bicgstab.f90 scalar_houc.f90 bicgstab2.f90
LIBFILE = libsolver.a

OFILES = $(F90FILES:.f90=.o)
MODFILES = $(F90FILES:.f90=.mod)

.SUFFIXES: .f90 .o

$(LIBFILE):$(OFILES) Makefile	
	cd $(OBJDIR); $(AR) $(LIBDIR)/$(LIBFILE) $(OFILES); $(RL) $(LIBDIR)/$(LIBFILE)

.f90.o: 
	$(F90) $(F90FLAGS) $(FLAGS) $(INCFLAGS) $(FFTW_INC) -c $*.f90 -o $(OBJDIR)/$*.o $(MODFLAGS)
 
clean: 
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(LIBDIR); rm -f $(LIBFILE)
