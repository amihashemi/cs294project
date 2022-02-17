#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "DBox.H"
#include "RectMDArray.H"
#include "FFT1D.H"
#include "FFTW1D.H"
#include "FFTMD.H"
#include "PoissonSolver.H"
#include "Projection.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
#include "RK4.H"
#include "WriteRectMDArray.H"
#include "PowerItoI.H"

int main(int argc, char* argv[])
{
  // Make N dependent on M (for safety)
  int M = 7;
  int N = Power(2, M);

  std::shared_ptr<FFT1D> fft1dptr(new FFTW1D(M));
  FieldData velocities(fft1dptr, 1);
  Projection testing(fft1dptr);

  int low[DIM] = {0,0};
  int high[DIM] = {N-1,N-1};
  DBox bx(low,high);
  double rho = 1.0/30.0; // careful when doing floating point arithmetic. don't want truncated int values
  double delta = 0.05;
  double dx = 1.0/((double)N); // careful when doing floating point arithmetic. don't want truncated int values

  for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt)) {
    if(pt[1]*dx <= 0.5) {
      velocities.m_data(pt, 0) = tanh((pt[1]*dx-0.25)/rho);
    } else if(pt[1]*dx > 0.5) {
      velocities.m_data(pt, 0) = tanh((0.75-pt[1]*dx)/rho);
    }

    velocities.m_data(pt, 1) = delta*sin(2.0*M_PI*pt[0]*dx);
  }

  testing.applyProjection(velocities);
  RectMDArray<double> divergenceSc(bx); // need to set the box for the MDArray
  testing.divergence(divergenceSc, velocities);

  MDWrite(string("ProjectionTest"), divergenceSc);

  ///somehow print out velocities
}

// Stuff from Projection:
//  Projection(std::shared_ptr<FFT1D> a_fft1dPtr);
//  void applyProjection(FieldData& a_velocity);

