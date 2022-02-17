// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#ifndef _POISSONSOLVER_H_
#define _POISSONSOLVER_H_

#include "FFT1D.H"
#include "FFTMD.H"
#include "RectMDArray.H"

class PoissonSolver
{
public:
  // Default constructor
  PoissonSolver();

  // Constructor
  PoissonSolver(std::shared_ptr< FFT1D > a_fft1dPtr);

  // Solve the Poisson term
  void solve(RectMDArray< double >& a_rhs) const;

  int m_M;
  int m_N;
  shared_ptr< FFT1D > m_fft1dPtr;
};

#endif

