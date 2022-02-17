//
// Created by Júlio Caineta on 15/12/17.
//

#include "DBox.H"
#include "FieldData.H"
#include "FFTW1D.H"
#include "WriteRectMDArray.H"

int main()
{
  int M = 7;
  int N = Power(2, M);
  int low[DIM] = {0,0};
  int high[DIM] = {N-1,N-1};
  DBox bx(low,high);

  std::shared_ptr<FFT1D> fft1dptr(new FFTW1D(M));
  FieldData velocities(fft1dptr, 0);
  velocities.m_data.setVal(9);
  velocities.setBoundaries();
  MDWrite(string("bc_test"), velocities.m_data);
};