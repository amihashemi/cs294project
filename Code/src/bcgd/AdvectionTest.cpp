#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Box.H"
#include "MDArray.H"
#include "FFT1D.H"
#include "FFTW1D.H"
#include "FFTMD.H"
#include "PoissonSolver.H"
#include "Projection.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
#include "RK4.H"
#include "WriteMDArray.H"  
#include "AdvectionOperator.H"
#include "PowerItoI.H"

int main(int argc, char* argv[])
{
  // Make N dependent on M (for safety)
  int M = 7;
  int N = Power(2, M);

  std::shared_ptr<FFT1D> fft1dptr(new FFTW1D(M));
  FieldData velocities(fft1dptr, 1);
  Projection testing(fft1dptr);
  DeltaVelocity divuu;
  divuu.init(velocities);
  int low[DIM] = {0,0};
  int high[DIM] = {N-1,N-1};
  Box bx(low,high);
  double rho = 1.0/30.0;
  double delta = 0.05;
  double dx = 1.0/((double)N);

  for(int i=0; i<bx.sizeOf(); i++) {
    int tuple[DIM];
    bx.tupleIndex(i, tuple);
      
      //printf("tuple 0 : %f tuple1 %f",tuple[0]*dx,tuple[1]*dx);

    if( (0.25 < tuple[0]*dx) && (tuple[0]*dx <= 0.5) && (0.25 < tuple[1]*dx) && (tuple[1]*dx <= 0.5)) {
        printf("tuple 0 : %f tuple1 %f",tuple[0]*dx,tuple[1]*dx);
        velocities.m_data[0][tuple] = 1.; //tanh((tuple[1]*dx-0.5)/rho);
         printf("vel : %f",velocities.m_data[0][tuple]);
    } else {
        velocities.m_data[0][tuple] = 0.5; //tanh((0.75-tuple[1]*dx)/rho);
    }

      velocities.m_data[1][tuple] = 0.5; //delta*sin(2.0*M_PI*tuple[0]*dx);
  }
    MDWrite("Advectioncheck[0]", velocities.m_data[0]);
    MDWrite("Advectioncheck[1]", velocities.m_data[1]);
    
  advectionOperator(divuu, velocities, bx, dx);

  MDWrite("Advection[0]", divuu.m_data[0]);
  MDWrite("Advection[1]", divuu.m_data[1]);

  return 0;
}
