#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "DBox.H"
#include "RectMDArray.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "FFTW1D.H"
#include "PoissonSolver.H"
#include "Projection.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
#include "RK4.H"
#include "WriteRectMDArray.H"
#include "ComputeEulerRHS.H"

double returnMax(RectMDArray<double, DIM> velocities) {
  double umax = 0.;
  double vmax = 0.;
  DBox grids = velocities[0].m_box;
  
  for(Point pt = grids.getLowCorner(); grids.notDone(pt); grids.increment(pt)) {
    double currVelocities[DIM];
    currVelocities[0] = velocities[0][pt];
    currVelocities[1] = velocities[1][pt];
    umax = max(currVelocities[0], umax);
    vmax = max(currVelocities[1], vmax);
  }

  return max(umax, vmax);
} 

void calcVorticity(RectMDArray<double, DIM> &vorticity, const DBox bx, FieldData &velocities, double dx) {
    velocities.fillGhosts();
  for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt)) {
    int tupleShiftLow[DIM];
    int tupleShiftHigh[DIM];
    double dudy;
    double dvdx;

    tupleShiftLow[0] = pt[0];
    tupleShiftLow[1] = pt[1]-1;
    tupleShiftHigh[0] = pt[0];
    tupleShiftHigh[1] = pt[1]+1;
    dudy = (velocities[0][tupleShiftHigh] - velocities[0][tupleShiftLow])/(2*dx);
    tupleShiftLow[0] = pt[0]-1;
    tupleShiftLow[1] = pt[1];
    tupleShiftHigh[0] = pt[0]+1;
    tupleShiftHigh[1] = pt[1];	
    dvdx = (velocities[1][tupleShiftHigh] - velocities[1][tupleShiftLow])/(2*dx);
    vorticity[tuple] = dudy-dvdx;
  }
}

int main(int argc, char* argv[])
{
  // Make N dependent on M (for safety)
  int M = 7;
  int N = Power(2, M);

  std::shared_ptr<FFT1D> fft1dptr(new FFTW1D(M));
  FieldData velocities(fft1dptr, 1);
  int low[DIM] = {0,0};
  int high[DIM] = {N-1,N-1};
  DBox bx(low,high);
  double rho = 1.0/30.0; // careful when doing floating point arithmetic. don't want truncated int values
  double delta = 0.05;
  double dx = 1.0/((double)N); // careful when doing floating point arithmetic. don't want truncated int values
  double timeEnd = 0.8;
  double C = 0.05;
  char name[10];
  int nstop = 15;

  /// Initial Conditions
  for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt)) {
    velocities.m_data[1][pt] = delta*sin(2.0*M_PI*pt[0]*dx);
    if(pt[1]*dx <= 0.5) {
      velocities.m_data[0][pt] = tanh((pt[1]*dx-0.25)/rho);
    } else if(pt[1]*dx > 0.5) {
      velocities.m_data[0][pt] = tanh((0.75-pt[1]*dx)/rho);
    }
  }

  RK4<FieldData, ComputeEulerRHS, DeltaVelocity> integrator;

  double dt = C*dx/returnMax(velocities.m_data);
  double time=0.0;

      RectMDArray<double, DIM> vorticity(bx);
      calcVorticity(vorticity, bx, velocities, dx);
      sprintf(name, "Vorticity.%.4f.%.4f", timeEnd, time);
      MDWrite(name, vorticity);
      int nstep = 0;
      while(nstep < nstop)
    {
      integrator.advance(time, dt, velocities);
      time += dt;
      nstep++;
      calcVorticity(vorticity, bx, velocities, dx);
      sprintf(name, "Vorticity.%d",nstep);
      MDWrite(name, vorticity);
      cout << "dt = " << dt << endl;
      dt = C*dx/returnMax(velocities.m_data);
    }
}
