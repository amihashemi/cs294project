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
  int M = 8;
  int N = Power(2, M);

  std::shared_ptr<FFT1D> fft1dptr(new FFTW1D(M));
  FieldData velocities(fft1dptr, 1);
  DeltaVelocity divuu;
  divuu.init(velocities);
  Point low = {{{0,0}}};
  Point high = {{{N-1,N-1}}};
  DBox bx(low, high);

  double dx = 1.0/((double)N);
    
    double sigma = 0.04;
    double x_0 = 0.4;
    double y_0 = 0.5;

    double sigma = 0.04;
    double x_0 = 0.4;
    double y_0 = 0.5;
    
    for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt)) {
        
        velocities.m_data(pt, 0) = exp(-(pow((pt[0]*dx-x_0),2)+
                                         pow((pt[1]*dx-y_0),2))/(2*pow(sigma,2.)))/pow(sigma,2.);
        
        velocities.m_data(pt, 1) = exp(-(pow((pt[0]*dx-x_0),2)+
                                         pow((pt[1]*dx-y_0),2))/(2*pow(sigma,2.)))/pow(sigma,2.);
    }
    
    
    MDWrite(string("velocitities"), velocities.m_data);

  advectionOperator(divuu, velocities, bx, dx);

  MDWrite(string("Advection"), divuu.m_data);

  return 0;
}
