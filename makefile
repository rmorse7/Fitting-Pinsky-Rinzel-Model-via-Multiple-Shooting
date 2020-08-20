CCC     = g++
f77	= f77
MAKRO  =  -DODESSA
#CFLAGS  = -g -O2  $(MAKRO)
#CCFLAGS = -g -O2  $(MAKRO)
CFLAGS  = -g  $(MAKRO)
CCFLAGS = -g  $(MAKRO)
INCDIR = -I/home/gabbiani/NAG/cll6i261dl/include
LIBS =   -lgfortran /usr/lib/x86_64-linux-gnu/libblas.so.3 ./libNAG/libnagc_nag.so ./libCVODES/libsundials_cvodes.a ./libCVODES/libsundials_nvecserial.a ./libCVODES/libsundials_shared.a ./libNAG/libifcoremt.so.5 ./libNAG/libimf.so ./libNAG/libirc.so ./libNAG/libsvml.so ./libNAG/libintlc.so.5 -lstdc++ -ldl -lpthread -lm

%.o: %.cc *.h makefile
	@echo -n [Compiling] $< 
	@$(CCC) $(CCFLAGS) $(INCDIR) -fPIC -o $@ -c $<
	@echo " ...done."

MODULES = modules/parse.o \
          modules/readData.o \
          modules/setMesh.o \
	  modules/initialise.o \
	  modules/simInit.o \
	  modules/setInitialValues.o \
	  modules/intODE.o \
	  modules/call_odessa.o \
	  modules/integrateRK.o \
	  modules/computeRight.o \
	  modules/solvLin.o \
	  modules/dampIt.o \
	  modules/outFit.o \
	  modules/freeMem.o

	  #modules/integrateCVODES.o CVODES not optimized!!

MODULES_SIMIT =	modules/parseSimit.o \
	        modules/readData.o \
		modules/setMesh.o \
		modules/initialise.o \
	        modules/simInit.o \
	  	modules/setInitialValues.o \
		modules/intODE.o \
		modules/call_odessa.o \
		modules/integrateRK.o \
		modules/computeRight.o \
		modules/solvLin.o \
		modules/dampIt.o \
		modules/outSimit.o \
		modules/freeMem.o

LIBSRES = 	libSRES/ESES.o \
		libSRES/ESSRSort.o \
		libSRES/sharefunc.o

MODEL = model.o

SUB = fitIt.o globOpt.o odessa.o dconstr.o setOptE04NCF.o lsei_wrapper.o

HEADERS = model.h def.h nr.h

NR = nr/nrutil.o \
     nr/rkck.o \
     nr/rkqs.o \
     nr/svdcmp.o \
     nr/svbksb.o \
     nr/pythag.o \
     nr/ran1.o \
     nr/gasdev.o \
     nr/spline.o \
     nr/choldc.o \
     nr/cholsl.o \
     nr/gam.o


all: dconstr.o odessa.o nr modules fitit diffit simit createSpline mex

diffit: diffit.o $(MODULES) $(LIBSRES)  $(SUB) $(MODEL) $(HEADERS) $(NR) makefile 
	@echo -n [Linking] $@
	@$(CCC) $(CCFLAGS) -o diffit diffit.o $(SUB) $(MODEL) $(NR) $(MODULES) $(LIBSRES) $(LIBS)
	@echo " ...done."


simit: simit.o $(MODULES_SIMIT) $(LIBSRES) $(SUB) $(MODEL) $(HEADERS) $(NR) makefile 
	@echo -n [Linking] $@
	@$(CCC) $(CCFLAGS) -o simit simit.o $(SUB) $(MODEL) $(NR) $(MODULES_SIMIT) $(LIBSRES) $(LIBS)
	@echo " ...done."

nr: $(NR)

# added -g flag for debugging
dconstr.o: dconstr.f makefile
	@echo -n [Compiling] $< 
	@$(f77) -g -fno-second-underscore -c dconstr.f -o dconstr.o
	@echo " ...done."

odessa.o: odessa.f makefile
	@echo -n [Compiling] $< 
	@$(f77) -fno-second-underscore -c odessa.f -o odessa.o
	@echo " ...done."

setOptE04NCF.o: setOptE04NCF.f makefile
	@echo -n [Compiling] $< 
	@$(f77) -fno-second-underscore -c setOptE04NCF.f -o setOptE04NCF.o
	@echo " ...done."

lsei_wrapper.o: lsei_wrapper.f90 makefile
	@echo -n [Compiling] $<
	@$(f77) -g -c lsei_wrapper.f90

modules: $(MODULES)

fitit: $(SUB)

model: $(MODEL)

createSpline:	createSpline.o $(NR) makefile 
	@echo -n [Linking] $@
	@$(CCC) $(CCFLAGS) -o createSpline createSpline.o $(NR)
	@echo " ...done."

mex:	
	@echo -n [Creating Mexfile] diffit_mex
	@mex -g diffit_mex.cc $(SUB) $(MODEL) $(NR) $(MODULES) $(LIBSRES) $(LIBS)
	@echo " ...done."
	@echo -n [Creating Mexfile] simit_mex
	@mex -g simit_mex.cc $(SUB) $(MODEL) $(NR) $(MODULES_SIMIT) $(LIBSRES) $(LIBS)
	@echo " ...done."
	@echo -n [Creating Mexfile] createSpline_mex
	@mex -g createSpline_mex.cc $(NR)
	@echo " ...done."


clean:
	@rm -f *.o
	@rm -f nr/*.o
	@rm -f modules/*.o
	@rm -f libSRES/*.o
	@rm -f *.mexglx

