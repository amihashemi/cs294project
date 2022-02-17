// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#include "RectMDArray.H"
#include "FFT1DW.H"
#include "FieldData.H"
#include "RK4.H"
#include "WriteRectMDArray.H"
#include "RHSNavierStokes.H"

void computeVorticity(RectMDArray< double >& vorticity, const DBox box,
                      FieldData& velocities, double dx)
{
  for (Point pt = box.getLowCorner(); box.notDone(pt); box.increment(pt))
  {
    Point lowShift;
    Point highShift;
    double dudy;
    double dvdx;

    lowShift[0] = pt[0];
    lowShift[1] = pt[1] - 1;
    highShift[0] = pt[0];
    highShift[1] = pt[1] + 1;
    dudy = (velocities.m_data(highShift, 0) -
            velocities.m_data(lowShift, 0)) / (2 * dx);
    lowShift[0] = pt[0] - 1;
    lowShift[1] = pt[1];
    highShift[0] = pt[0] + 1;
    highShift[1] = pt[1];
    dvdx = (velocities.m_data(highShift, 1) -
            velocities.m_data(lowShift, 1)) / (2 * dx);
    vorticity[pt] = dudy - dvdx;
  }
};

double maxVelocity(RectMDArray< double, DIM > velocities)
{
  // only working for 2D
  double umax = 0.;
  double vmax = 0.;
  DBox box = velocities.getDBox();

  for (Point pt = box.getLowCorner(); box.notDone(pt); box.increment(pt))
  {
    double velocity[2];
    velocity[0] = velocities(pt, 0);
    velocity[1] = velocities(pt, 1);
    umax = max(velocity[0], umax);
    vmax = max(velocity[1], vmax);
  }
  return max(umax, vmax);
};

int main(int argc, char* argv[])
{
  int M;
//  int M = 8;
  cout << "input log_2(number of grid points) [e.g., 8]" << endl;
  cin >> M;
  int N = Power(2, M);

  std::shared_ptr< FFT1D > fft1dPtr(new FFT1DW(M));
  FieldData velocities(fft1dPtr, 1);
  Point low = {{{0, 0}}};
  Point high = {{{N - 1, N - 1}}};
  DBox box(low, high);

  double dx = 1.0 / ((double) N);
  double timeEnd = 0.8;
  char filename[10];
  int nstop = 400;

  double sigma = 0.04;
  double x_0 = 0.4;
  double y_0 = 0.5;

  for (Point pt = box.getLowCorner(); box.notDone(pt); box.increment(pt))
  {

    velocities.m_data(pt, 0) = exp(-(pow((pt[0] * dx - x_0), 2) +
                                     pow((pt[1] * dx - y_0), 2)) /
                                   (2 * pow(sigma, 2.))) / pow(sigma, 2.);

    velocities.m_data(pt, 1) = exp(-(pow((pt[0] * dx - x_0), 2) +
                                     pow((pt[1] * dx - y_0), 2)) /
                                   (2 * pow(sigma, 2.))) / pow(sigma, 2.);
  }

  RK4 <FieldData, RHSNavierStokes, DeltaVelocity> integrator;

  double dt = 0.0001;
  double time = 0.0;

  RectMDArray< double > vorticity(box);
  computeVorticity(vorticity, box, velocities, dx);
  sprintf(filename, "Vorticity.%.4f.%.4f", timeEnd, time);
  MDWrite(string(filename), vorticity);
  int nstep = 0;
  while (nstep < nstop)
  {
    velocities.setBoundaries();
    integrator.advance(time, dt, velocities);
    sprintf(filename, "velocity.%d", nstep);
    MDWrite(string(filename), velocities.m_data);
    velocities.setBoundaries();
    time += dt;
    nstep++;
    computeVorticity(vorticity, box, velocities, dx);
    sprintf(filename, "Vorticity.%d", nstep);
    MDWrite(string(filename), vorticity);
    cout << "iter = " << nstep << endl;
  }
};
