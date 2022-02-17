// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#include "DBox.H"
#include "FieldData.H"
#include "FFT1DW.H"
#include "WriteRectMDArray.H"

int main()
{
  int M = 7;
  int N = Power(2, M);
  int low[2] = {0, 0};
  int high[2] = {N - 1, N - 1};
  DBox bx(low, high);

  std::shared_ptr< FFT1D > fft1dptr(new FFT1DW(M));
  FieldData velocities(fft1dptr, 0);
  velocities.m_data.setVal(9);
  velocities.setBoundaries();
  MDWrite(string("bc_test"), velocities.m_data);
};