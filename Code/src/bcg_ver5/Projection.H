#ifndef _PROJECTION_H_
#define _PROJECTION_H_

#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "DBox.H"
#include "RectMDArray.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
#include "PoissonSolver.H"

class Projection
{
public:
  //Projection();
  Projection(std::shared_ptr<FFT1D> a_fft1dPtr);
  //~Projection();
  void applyProjection(FieldData& a_velocity);
  void gradient(DeltaVelocity& a_vector,const RectMDArray<double>& a_scalar) const;
  void divergence(RectMDArray<double>& a_scalar,const FieldData& a_vector) const;

  int m_N;
  //int m_M;    // these turned out to be unnecessary
  //DBox m_grid; // these turned out to be unnecessary
  PoissonSolver m_solver;
};

#endif
