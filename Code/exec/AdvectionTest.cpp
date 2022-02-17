// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#include "RectMDArray.H"
#include "FFT1DW.H"
#include "PoissonSolver.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
#include "WriteRectMDArray.H"
#include "AdvectionSolver.H"
#include "PowerItoI.H"

int main(int argc, char* argv[])
{
  int M;
//  int M = 8;
  cout << "input log_2(number of grid points) [e.g., 8]" << endl;
  cin >> M;
  int N = Power(2, M);

  std::shared_ptr< FFT1D > fft1dPtr(new FFT1DW(M));
  FieldData velocities(fft1dPtr, 1);
  DeltaVelocity divuu;
  divuu.init(velocities);
  Point low = {{{0, 0}}};
  Point high = {{{N - 1, N - 1}}};
  DBox bx(low, high);

  double dx = 1.0 / ((double) N);

  double sigma = 0.04;
  double x_0 = 0.4;
  double y_0 = 0.5;

  for (Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt))
  {

    velocities.m_data(pt, 0) = exp(-(pow((pt[0] * dx - x_0), 2) +
                                     pow((pt[1] * dx - y_0), 2)) /
                                   (2 * pow(sigma, 2.))) / pow(sigma, 2.);

    velocities.m_data(pt, 1) = exp(-(pow((pt[0] * dx - x_0), 2) +
                                     pow((pt[1] * dx - y_0), 2)) /
                                   (2 * pow(sigma, 2.))) / pow(sigma, 2.);
  }

  advectionSolver(divuu, velocities, bx, dx);
  MDWrite(string("velocitities"), velocities.m_data);
  MDWrite(string("Advection"), divuu.m_data);

  return 0;
};
