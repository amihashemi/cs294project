#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <list>
#include <memory>
#include "DBox.H"
#include "RectMDArray.H"
#include "FFT1D.H"
#include "FFTMD.H"
using namespace std;
FFTMD::FFTMD(std::shared_ptr<FFT1D> a_fft1dPtr)
{
  m_fft1dPtr = a_fft1dPtr;
  m_M = m_fft1dPtr->getM();
  m_N = m_fft1dPtr->getN();
} 
void FFTMD::forwardCC(RectMDArray<complex<double> > & a_f) const
{
//  int low[DIM],high[DIM],tuple[DIM];
  int low[DIM],high[DIM];
  vector<complex<double> > f1d(m_N);
  vector<complex<double> > fHat1d(m_N);
  
  for (int dir = 0;dir < DIM ; dir++)
    {
      for (int dir2 = 0;dir2 < DIM;dir2++)
        {
          low[dir2]= 0;
          high[dir2] = m_N-1;
        }
      high[dir]=0;
      DBox base(low,high);
//      for (int k = 0;k < base.sizeOf();k++)
      for(Point pt = base.getLowCorner(); base.notDone(pt); base.increment(pt))
        {
//          base.tupleIndex(k,tuple);
          for (int l = 0 ; l < m_N;l++)
            {
              pt[dir]= l;
              f1d[l] = a_f[pt];
            }
          m_fft1dPtr->forwardFFTCC(fHat1d,f1d);
          pt[dir] = 0;
          for (int l = 0 ; l < m_N;l++)
            {
              pt[dir] = l;
              a_f[pt] = fHat1d[l];
            }
        }
    }     
};
void FFTMD::inverseCC(RectMDArray<complex<double> > & a_fHat) const
{int low[DIM],high[DIM];
  vector<complex<double> > f1d(m_N);
  vector<complex<double> > fHat1d(m_N);
  
  for (int dir = 0;dir < DIM ; dir++)
    {
      for (int dir2 = 0;dir2 < DIM;dir2++)
        {
          low[dir2]= 0;
          high[dir2] = m_N-1;
        }
      high[dir]=0;
      DBox base(low,high);
      for(Point pt = base.getLowCorner(); base.notDone(pt); base.increment(pt))
//      for (int k = 0;k < base.sizeOf();k++)
        {
//          base.tupleIndex(k,tuple);
          for (int l = 0 ; l < m_N;l++)
            {
              pt[dir] = l;
              fHat1d[l] = a_fHat[pt];
            }
          m_fft1dPtr->inverseFFTCC(f1d,fHat1d);
          pt[dir] = 0;
          for (int l = 0 ; l < m_N;l++)
            {
              pt[dir] = l;
              a_fHat[pt] = f1d[l];
            }
        }
    }     
};
const int& FFTMD::getN() const
{
  return m_N;
  
};
const int& FFTMD::getM() const
{
  return m_M;
};
