include ../Makefile.in

F90FILES = precision.f90 string.f90 cli_reader.f90 memory.f90 fileio.f90 parser.f90 infnan.f90 parallel.f90 \
	   random.f90 math.f90 pentadiagonal.f90 polydiagonal.f90 tridiagonal.f90 quicksort.f90 quadrature.f90 \
	utilities.f90
CFILES = fileio_c.c
LIBFILE = liblibrary.a

OFILES = $(F90FILES:.f90=.o) $(CFILES:.c=.o)
MODFILES = $(F90FILES:.f90=.mod)

.SUFFIXES: .o .f90 .c

$(LIBFILE):$(OFILES) Makefile	
	cd $(OBJDIR); $(AR) $(LIBDIR)/$(LIBFILE) $(OFILES); $(RL) $(LIBDIR)/$(LIBFILE)

.f90.o: 
	$(F90) $(F90FLAGS) $(FLAGS) $(INCFLAGS) -c $*.f90 -o $(OBJDIR)/$*.o $(MODFLAGS)
.c.o: 
	$(CC) $(CFLAGS) $(INCFLAGS) -c $*.c -o $(OBJDIR)/$*.o
 
clean: 
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(LIBDIR); rm -f $(LIBFILE)
