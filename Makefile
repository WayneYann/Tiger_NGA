include Makefile.in

# Files
CODE = arts
MODULES = library geometry metric core monitor operator solver combustion sgsmodel boundary pollutants spray radiation postprocess
LIBRARIES = core postprocess sgsmodel combustion solver radiation spray pollutants boundary operator geometry metric monitor library

# Targets

default:Makefile Makefile.in
	@make libraries
	@make $(CODE)
	@make util

all:    Makefile Makefile.in 
	@make libraries
	@make $(CODE)
	@make util

opt:
	@make "FLAGS = $(OPTFLAGS)"

debug:
	@make "FLAGS = $(DBGFLAGS)"

allopt:
	@make all "FLAGS = $(OPTFLAGS)"

alldebug:
	@make all "FLAGS = $(DBGFLAGS)"

$(CODE):$(MODULES:%=lib%.a)
	$(LD) $(FLAGS) $(OBJDIR)/main.o $(LIBRARIES:%=$(LIBDIR)/lib%.a) $(LAPACK_LIB) $(BLAS_LIB) $(HYPRE_LIB) $(FFTW_LIB) $(SUNDIALS_LIB) -o $(BINDIR)/$(CODE) $(LDFLAGS)

libraries:
	@for i in $(MODULES); do make -C $$i; done

clean:
	@for i in $(MODULES); do make clean -C $$i; done
	rm -f $(BINDIR)/$(CODE)

util:  $(MODULES:%=lib%.a)
	@make -C tools

utilclean:
	@make -C tools clean

distclean:
	@make clean
	@make utilclean

install:
	$(INSTSCPT)
