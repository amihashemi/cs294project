#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

system := $(shell uname)
FFTWDIR = /usr/local
FFTINCFLAGS = -I$(FFTWDIR)/include
#for Johnny and possibly Ivan (Ubuntu machines)
ifeq ((system),Linux)
  #LAPACKFLAGS = -lcblas -lclapack -I/usr/include/atlas -L/usr/local/lib -L/usr/lib
  LAPACKFLAGS = -lcblas -llapack -I/usr/include/atlas
endif
# for OSX
ifeq ($(system),Darwin)
  LAPACKFLAGS = -framework Accelerate
endif
# don't actually need LAPACK
LAPACKFLAGS = 

DIM = 2
FFTTOOLS = ../fftTools
RMDA = ../RectMDArray
TIMER = ../../utils/timer
VPATH = ../../utils/Writers
CFLAGS = -g -Wall $(LAPACKFLAGS)
CFLAGS += -std=c++11
CFLAGS += -DDIM=$(DIM)
#CFLAGS = -O3
LIBS:=  -L$(FFTWDIR)/lib -lfftw3 
#CXX = g++
CXX = clang++

CFLAGS += $(FFTINCFLAGS) -I$(VPATH) -I$(RMDA) -I$(FFTTOOLS) -I$(TIMER)

REQOBJS= FFTMD.o FFTW1D.o PowerItoI.o DBox.o $(VPATH)/WriteRectMDArray.o PoissonSolver.o VisitWriter.o

EULEROBJS= Projection.o AdvectionOperator.o ComputeEulerRHS.o

ProjectionTestObjs = $(REQOBJS) FieldData.o DeltaVelocity.o Projection.o

AdvectionTestObjs = $(REQOBJS) FieldData.o DeltaVelocity.o AdvectionOperator.o

ComputeEulerRHSTestObjs = $(REQOBJS) FieldData.o DeltaVelocity.o Projection.o AdvectionOperator.o ComputeEulerRHS.o

BCTestObjs = $(REQOBJS) FieldData.o

FFTMD.o: FFTMD.cpp FFTMD.H $(RMDA)/DBox.H $(RMDA)/RectMDArray.H $(RMDA)/RectMDArrayImplem.H GNUmakefile
	$(CXX) $(CFLAGS) -c FFTMD.cpp

FFTW1D.o: FFTW1D.H FFTW1D.cpp GNUmakefile
	$(CXX) $(CFLAGS) $(FFTINCFLAGS) -c -o FFTW1D.o FFTW1D.cpp

PowerItoI.o: PowerItoI.cpp PowerItoI.H GNUmakefile
	$(CXX) $(CFLAGS) -c PowerItoI.cpp

DBox.o: $(RMDA)/DBox.cpp $(RMDA)/DBox.H GNUmakefile
	$(CXX) $(CFLAGS) -c $(RMDA)/DBox.cpp

WriteRectMDArray.o: $(VPATH)/WriteRectMDArray.cpp $(VPATH)/WriteRectMDArray.H $(RMDA)/RectMDArray.H $(RMDA)/RectMDArrayImplem.H VisitWriter.H GNUmakefile
	$(CXX) $(CFLAGS) -c $(VPATH)/WriteRectMDArray.cpp

DeltaVelocity.o: DeltaVelocity.H DeltaVelocity.cpp FieldData.H $(RMDA)/DBox.H $(RMDA)/RectMDArray.H $(RMDA)/RectMDArrayImplem.H GNUmakefile
	$(CXX) $(CFLAGS) -c DeltaVelocity.cpp

FieldData.o: FieldData.H FieldData.cpp DeltaVelocity.H $(RMDA)/DBox.H $(RMDA)/RectMDArray.H $(RMDA)/RectMDArrayImplem.H FFT1D.H GNUmakefile
	$(CXX) $(CFLAGS) -c FieldData.cpp

PoissonSolver.o: PoissonSolver.H PoissonSolver.cpp $(RMDA)/DBox.H FFT1D.H FFTMD.H GNUmakefile
	$(CXX) $(CFLAGS) -c PoissonSolver.cpp

Projection.o: Projection.H Projection.cpp PoissonSolver.H $(RMDA)/DBox.H $(RMDA)/RectMDArray.H $(RMDA)/RectMDArrayImplem.H FieldData.H DeltaVelocity.H GNUmakefile
	$(CXX) $(CFLAGS) -c Projection.cpp

AdvectionOperator.o: AdvectionOperator.H AdvectionOperator.cpp $(RMDA)/DBox.H $(RMDA)/RectMDArray.H FieldData.H DeltaVelocity.H GNUmakefile
	$(CXX) $(CFLAGS) -c AdvectionOperator.cpp

MyMachineName:
	@echo $(system) $(HOME)

AdvectionTest: $(AdvectionTestObjs) AdvectionTest.cpp GNUmakefile
	$(CXX) $(CFLAGS) AdvectionTest.cpp $(AdvectionTestObjs) -o AdvectionTest $(LIBS)

ProjectionTest: $(ProjectionTestObjs) ProjectionTest.cpp GNUmakefile
	$(CXX) $(CFLAGS) ProjectionTest.cpp $(ProjectionTestObjs) -o ProjectionTest $(LIBS)

EulerTest: $(ComputeEulerRHSTestObjs) EulerTest.cpp $(FFTTOOLS)/RK4.H GNUmakefile
	$(CXX) $(CFLAGS) EulerTest.cpp $(ComputeEulerRHSTestObjs) -o EulerTest $(LIBS)

BCTest: $(BCTestObjs) bc_test.cpp GNUmakefile
	$(CXX) $(CFLAGS) bc_test.cpp $(BCTestObjs) -o bc_test $(LIBS)

%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

clean:
	rm *.o *.d
cleanExecs:
	rm AdvectionTest ProjectionTest EulerTest
