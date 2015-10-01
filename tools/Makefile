include ../Makefile.in

BINFILES = init_flow data2plot3d config2plot3d chucks2arts \
	   config2ensight data2ensight createChemtable \
	   dns2les interpolatedata \
	   biinflow2inflow chemtable2plot3d \
	   editData statdata2plot3d statconfig2plot3d spectrum \
	   rescaler compare mergeData bl2ml channel2jet \
	   create1DChemtable combineInflow createPremtable \
           createUnsteadyPremtable

all:	Makefile
	@for i in $(BINFILES); do make -C $$i; done

clean: 
	@for i in $(BINFILES); do make clean -C $$i; done

# premtable2plot3d