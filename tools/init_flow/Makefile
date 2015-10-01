include ../../Makefile.in

F90FILES = param.f90 taylor.f90 init_flow.f90 dns_box.f90 channel.f90 \
	   pipe.f90 jetcyl.f90 jetcart.f90 mixing.f90 gambit.f90 \
	   gambit_reader.f90 conv_scalar.f90 conv_vortex.f90 \
	   spatial_mixing.f90 laminar_flame.f90 vonkarman.f90 \
	   lamb_vortex.f90 density_vortex.f90 \
	   boundary_layer.f90 \
	   gfm_test.f90 marmottant.f90 \
	   milk_crown.f90 deformation.f90 \
	   rayleigh_taylor.f90 round_cart_jet.f90 diesel_dns.f90 \
           energy_cons.f90 vortex_ring.f90 \
	   rt_vortex.f90 \
	   rayleigh.f90 wave.f90 spurious_currents.f90 curvature.f90 \
	   ignition.f90 blowing.f90 air_layer.f90 \
	   helium_plume.f90 grating.f90 bubble.f90 \
           drop_coll.f90 glass.f90 density_ratio.f90 uniform.f90 \
	   jetcart2.f90 simplejet.f90 jetcyl_coflow.f90


BINFILE = init_flow

OFILES   = $(F90FILES:.f90=.o)
MODFILES = $(F90FILES:.f90=.mod)

.SUFFIXES: .f90 .o

$(BINFILE): $(OFILES) $(LIBDIR)/liblibrary.a $(LIBDIR)/libcombustion.a
	cd $(OBJDIR);$(LD) $(FLAGS) -o $(BINDIR)/$(BINFILE) $+ \
	$(FFTW_LIB) $(LAPACK_LIB) $(BLAS_LIB) $(LDFLAGS) 

clean: 
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(BINDIR); rm -f $(BINFILE)

.f90.o: 
	$(F90) $(F90FLAGS) $(FLAGS) $(INCFLAGS) $(FFTW_INC) -c $*.f90 -o $(OBJDIR)/$*.o $(MODFLAGS)
