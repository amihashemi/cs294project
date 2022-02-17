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
  DBox grids = velocities.getDBox();
  
  for(Point pt = grids.getLowCorner(); grids.notDone(pt); grids.increment(pt)) {
    double currVelocities[DIM];
    currVelocities[0] = velocities(pt, 0);
    currVelocities[1] = velocities(pt, 1);
    umax = max(currVelocities[0], umax);
    vmax = max(currVelocities[1], vmax);
  }

  return max(umax, vmax);
} 

void calcVorticity(RectMDArray<double> &vorticity, const DBox bx, FieldData &velocities, double dx) {
//    velocities.fillGhosts();
  for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt)) {
    Point tupleShiftLow;
    Point tupleShiftHigh;
    double dudy;
    double dvdx;

    tupleShiftLow[0] = pt[0];
    tupleShiftLow[1] = pt[1]-1;
    tupleShiftHigh[0] = pt[0];
    tupleShiftHigh[1] = pt[1]+1;
    dudy = (velocities.m_data(tupleShiftHigh, 0) - velocities.m_data(tupleShiftLow, 0))/(2*dx);
    tupleShiftLow[0] = pt[0]-1;
    tupleShiftLow[1] = pt[1];
    tupleShiftHigh[0] = pt[0]+1;
    tupleShiftHigh[1] = pt[1];	
    dvdx = (velocities.m_data(tupleShiftHigh, 1) - velocities.m_data(tupleShiftLow, 1))/(2*dx);
    vorticity[pt] = dudy-dvdx;
//    vorticity(pt, 0) = dudy-dvdx; // this is a workaround: do we really need 2D here?
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

    double sigma = 0.02;
    double x_0 = 0.4;
    double y_0 = 0.5;
    
    //  for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt)) {
    //    if(pt[1]*dx <= 0.5) {
    //      velocities.m_data(pt, 0) = tanh((pt[1]*dx-0.25)/rho);
    //    } else if(pt[1]*dx > 0.5) {
    //      velocities.m_data(pt, 0) = tanh((0.75-pt[1]*dx)/rho);
    //    }
    //
    //    velocities.m_data(pt, 1) = delta*sin(2.0*M_PI*pt[0]*dx);
    //  }
    //
    for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt)) {
        //        if(pt[1]*dx <= 0.5) {
        velocities.m_data(pt, 0) = exp(-(pow((pt[0]*dx-x_0),2)+
                                         pow((pt[1]*dx-y_0),2))/(2*pow(sigma,2.)))/pow(sigma,2.);
        //        } else if(pt[1]*dx > 0.5) {
        //            velocities.m_data(pt, 0) = tanh((0.75-pt[1]*dx)/rho);
        //        }
        
        velocities.m_data(pt, 1) = exp(-(pow((pt[0]*dx-x_0),2)+
                                         pow((pt[1]*dx-y_0),2))/(2*pow(sigma,2.)))/pow(sigma,2.);
    }
    
  RK4<FieldData, ComputeEulerRHS, DeltaVelocity> integrator;

    double dt = 0.0000001; //C*dx/returnMax(velocities.m_data);
  double time=0.0;

      string outnames[2] = {{"Vorticity[0]", "Vorticity[1]"}};
      array<double, DIM> corner;
      corner.fill(0.);
      double h = 1./(bx.getHighCorner()[0]);

      RectMDArray<double> vorticity(bx);
      calcVorticity(vorticity, bx, velocities, dx);
      sprintf(name, "Vorticity.%.4f.%.4f", timeEnd, time);
      MDWrite(name, vorticity, outnames, corner, h);
      int nstep = 0;
      while(nstep < nstop)
    {
        velocities.setBoundaries();
      integrator.advance(time, dt, velocities);
        velocities.setBoundaries();
      time += dt;
      nstep++;
      calcVorticity(vorticity, bx, velocities, dx);
        sprintf(name, "velocity.%d",nstep);
          MDWrite(name, velocities.m_data, outnames, corner, h);
      sprintf(name, "Vorticity.%d",nstep);
      MDWrite(name, vorticity, outnames, corner, h);
      cout << "dt = " << dt << endl;
        dt = 0.0000001;  //C*dx/returnMax(velocities.m_data);
    }
}
