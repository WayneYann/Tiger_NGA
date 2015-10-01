include ../../Makefile.in

#F90FILES = flamelet_4D.f90 table_4D.f90 createChemtable_4D.f90
F90FILES = unsteady_flamelet_tools.f90 unsteady_flamelet.f90 unsteady_table.f90 createUnsteadyChemtable.f90

BINFILE = createUnsteadyChemtable

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
