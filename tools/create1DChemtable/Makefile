include ../../Makefile.in

F90FILES = flamelet_1D.f90 table_1D.f90 create1DChemtable.f90

BINFILE = create1DChemtable

OFILES   = $(F90FILES:.f90=.o)
MODFILES = $(F90FILES:.f90=.mod)

.SUFFIXES: .f90 .o

$(BINFILE): $(OFILES) $(LIBDIR)/liblibrary.a
	cd $(OBJDIR);$(LD) $(FLAGS) -o $(BINDIR)/$(BINFILE) $+ $(LDFLAGS)

clean: 
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(BINDIR); rm -f $(BINFILE)

.f90.o: 
	$(F90) $(F90FLAGS) $(FLAGS) $(INCFLAGS) -c $*.f90 -o $(OBJDIR)/$*.o $(MODFLAGS)
