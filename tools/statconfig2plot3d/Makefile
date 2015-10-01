include ../../Makefile.in

F90FILES = statconfig2plot3d.f90

BINFILE = statconfig2plot3d

OFILES   = $(F90FILES:.f90=.o)
MODFILES = $(F90FILES:.f90=.mod)

.SUFFIXES: .f90 .o

$(BINFILE): $(LIBDIR)/liblibrary.a $(OFILES)
	cd $(OBJDIR);$(LD) $(LDFLAGS) $(FLAGS) $(INCFLAGS) $(LIBFLAGS) -o $(BINDIR)/$(BINFILE) $+ $(LIBDIR)/liblibrary.a

clean: 
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(BINDIR); rm -f $(BINFILE)

.f90.o: 
	$(F90) $(F90FLAGS) $(FLAGS) $(INCFLAGS) $(FFTW_INC) -c $*.f90 -o $(OBJDIR)/$*.o $(MODFLAGS)
