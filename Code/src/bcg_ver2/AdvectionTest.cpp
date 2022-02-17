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
//#include "Projection.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
#include "RK4.H"
#include "WriteRectMDArray.H"
#include "AdvectionOperator.H"
#include "PowerItoI.H"

int main(int argc, char* argv[])
{
  // Make N dependent on M (for safety)
  int M = 7;
  int N = Power(2, M);

  std::shared_ptr<FFT1D> fft1dptr(new FFTW1D(M));
  FieldData velocities(fft1dptr, 1);
//  Projection testing(fft1dptr);
  DeltaVelocity divuu;
  divuu.init(velocities);
  int lowpt[DIM] = {0,0};
  int highpt[DIM] = {N-1,N-1};
  Point low(lowpt);
  Point high(highpt);
  DBox bx(low, high);
  double rho = 1.0/30.0;
  double delta = 0.05;
  double dx = 1.0/((double)N);

  for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt)) {
    if(pt[1]*dx <= 0.5) {
        velocities.m_data(pt, 1) = tanh((pt[1]*dx-0.25)/rho);
    } else if(pt[1]*dx > 0.5) {
        velocities.m_data(pt, 1) = tanh((0.75-pt[1]*dx)/rho);
    }

    velocities.m_data(pt, 0) = delta*sin(2.0*M_PI*pt[0]*dx);
  }
 
  advectionOperator(divuu, velocities, bx, dx);
 //a_divuu.m_data(Point(accessTuple), int(!k)) = a_divuu.m_data(Point(accessTuple), k) + deriv;
//  MDWrite("Advection[0]", divuu.m_data[0]);
//  MDWrite("Advection[1]", divuu.m_data[1]);
  string outnames[2] = {{"Advection[0]", "Advection[1]"}};
  array<double, DIM> corner;
  corner.fill(0.);
  double h = 1./(highpt[0]);
    MDWrite("Velocity", velocities.m_data, outnames, corner, h);
  MDWrite("Advection", divuu.m_data, outnames, corner, h);

  return 0;
}
