#ifndef _POISSONSOLVER_H_
#define _POISSONSOLVER_H_

#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <complex>
#include "DBox.H"
#include "RectMDArray.H"
#include "FFT1D.H"
#include "FFTMD.H"
using namespace std;

class PoissonSolver
{
public:
  PoissonSolver();
  PoissonSolver(std::shared_ptr<FFT1D> a_fft1dPtr);
  void solve(RectMDArray<double>& a_Rhs) const;
  int m_M,m_N;
  std::shared_ptr<FFT1D> m_fft1dptr;
};

#endif

