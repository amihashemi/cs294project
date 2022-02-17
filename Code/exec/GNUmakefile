#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
HOME = ./..
WRITERS = $(HOME)/utils/Writers
RMDA = $(HOME)/src/RectMDArray
FFT = $(HOME)/src/fftTools
BCHLAI = $(HOME)/src/bchlai
# BCHLAI = $(HOME)/src/bcg
TIMER = $(HOME)/utils/timer
VPATH = . $(HOME) $(BCHLAI) $(RMDA) $(FFT) $(TIMER) $(WRITERS)

# set problem dimensions
DIM = 2

#CXX = g++
CXX = clang++
FFTWDIR = /usr/local
CFLAGS = -g -Wall
#CFLAGS = -O3
CFLAGS += -std=c++11 -I$(BCHLAI) -I$(RMDA) -I$(FFT) -I$(TIMER) -I$(WRITERS)
CFLAGS += -DDIM=$(DIM)
LIBS:= -L$(FFTWDIR)/lib -lfftw3

OBJS = $(FFT)/FFTMD.o $(FFT)/FFT1DW.o $(FFT)/PowerItoI.o $(RMDA)/DBox.o $(WRITERS)/WriteRectMDArray.o $(WRITERS)/VisitWriter.o

AdvectionTestObjs = $(OBJS) AdvectionTest.cpp $(BCHLAI)/FieldData.o $(BCHLAI)/DeltaVelocity.o $(BCHLAI)/AdvectionSolver.o

BCHLaiTestObjs = $(OBJS) BCHLaiTest.cpp $(BCHLAI)/FieldData.o $(BCHLAI)/DeltaVelocity.o $(BCHLAI)/Projection.o $(BCHLAI)/AdvectionSolver.o $(BCHLAI)/RHSNavierStokes.o $(BCHLAI)/PoissonSolver.o

BCTestObjs = $(OBJS) BCTest.cpp $(BCHLAI)/FieldData.o


%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

BCHLaiTest: GNUmakefile $(BCHLaiTestObjs)  $(FFT)/RK4.H
	$(CXX) $(CFLAGS) $(BCHLaiTestObjs) $(LIBS) -o BCHLaiTest$(DIM)D.exe

AdvectionTest: GNUmakefile $(AdvectionTestObjs)
	$(CXX) $(CFLAGS) $(AdvectionTestObjs) $(LIBS) -o AdvectionTest$(DIM)D.exe

BCTest: GNUmakefile $(BCTestObjs)
	$(CXX) $(CFLAGS) $(BCTestObjs) $(LIBS) -o bCTest.exe

clean:
	rm -r *.exe *.o *.d

clean_vtk:
	rm *.vtk
