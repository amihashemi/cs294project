#include <iostream>
#include "PoissonSolver.H"
using namespace std;

PoissonSolver::PoissonSolver()
{
}
PoissonSolver::PoissonSolver(std::shared_ptr<FFT1D> a_fft1dPtr)
{
  m_fft1dptr = a_fft1dPtr;
  m_N = a_fft1dPtr->getN();
  m_M = a_fft1dPtr->getM();
}

void PoissonSolver::solve(RectMDArray<double>& a_Rhs) const {
  //create a contructor
  DBox bx = a_Rhs.getDBox();
  RectMDArray<complex<double> > fftwForward(bx);
  RectMDArray<complex<double> > fourierCoef(bx);

  double h = 1.0/((double)m_N);
   h = .0001;
  //double h = 1/m_N; // integer division, becomes 0.0 in the end!

//  for (int k = 0; k < bx.sizeOf();k++) {
  for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt)) {
    fftwForward[pt] = complex<double>(a_Rhs[pt], 0.0);
  }
  FFTMD fftmd = FFTMD(m_fft1dptr);
  fftmd.forwardCC(fftwForward);

//  for (int k=0; k<bx.sizeOf(); k++) {
  for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt)) {
//    int index[2];
//    bx.tupleIndex(k,index);
    int i = pt[0];
    int j = pt[1];
    if (pt == bx.getLowCorner()) {
      fourierCoef[pt] = complex<double>(0.0,0.0);
    } else {
      complex<double> div((2.0*cos(2.0*M_PI*((double)i)*h) - 2.0*cos(2.0*M_PI*((double)j)*h) - 4.0)/h/h,0.0);
      fourierCoef[pt] =fftwForward[pt]/div;
    }
  }

  fftmd.inverseCC(fourierCoef);

//  for (int k=0; k<bx.sizeOf(); k++) {
  for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt)) {
    a_Rhs[pt] = real(fourierCoef[pt]) * h * h;
  }

}
