#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "DBox.H"
#include "RectMDArray.H"
#include "FFTW1DSine.H"
#include "FFTMDSine.H"
#include "PoissonSolver.H"
#include "Projection.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
#include "RK4.H"
#include "WriteRectMDArray.H"

int main(int argc, char* argv[])
{
  unsigned int N = 128;
  int M = log(N)/log(2);
  std::shared_ptr<FFT1D> fft1dptr(M);
  FieldData velocities(fft1dptr, 1);
  PoissonSolver testing(fft1dptr);
  int low[DIM] = {0,0};
  int high[DIM] = {N,N};
  DBox bx(low,high);
  double rho = 1/30;
  double delta = 0.05;
  double dx = 1/N;
  double timeEnd = 0.8;
  double C = 0.5;
  char name[10];

  velocities.m_grid = bx;
  velocities.m_N = N;
  velocities.m_M = log(N)/log(2);
  for(int i=0; i<DIM; i++) {
    velocities[i].define(bx);
  }

  for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt)) {
    if(pt[1]*dx <= 0.5) {
      velocities.m_data[pt][0] = tanh((y-0.25)/rho);
    } else if(pt[1]*dx > 0.5) {
      velocities.m_data[pt][0] = tanh((0.75-y)/rho);
    }

    velocities.m_data[pt][1] = delta*sin(2*M_PI*pt[0]*dx);
  }

//  testing.solve(m_data[0]);
//  testing.solve(m_data[1]);
  testing.solve(m_data);

  ///somehow print out velocities
}
