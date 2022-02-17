// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#ifndef _PROJECTION_H
#define _PROJECTION_H

#include "PoissonSolver.H"
#include "RectMDArray.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "FieldData.H"
#include "DeltaVelocity.H"

/// Implements the projection approximation according to Lai (1993)
class Projection
{
public:
  // Default constructor
  Projection();

  // Constructor
  Projection(int a_M);

  Projection(std::shared_ptr< FFT1D > a_fft1dPtr);

  // Default destructor
  ~Projection(){};

  // Apply the projection operator to a given vector field (velocity)
  void applyProjection(FieldData& a_velocity) const;

  // Compute the gradient of a given scalar field
  void gradient(DeltaVelocity& a_vector,
                const RectMDArray< double >& a_scalar) const;

  // Compute the divergence of a given vector field
  void divergence(RectMDArray< double >& a_scalar,
                  const FieldData& a_vector) const;

private:
  int m_M;
  int m_N;
  DBox m_grid;
/**
  The divergence of the velocity field can be written as a Poisson equation,
  which is solved for the scalar Phi in the first step of the computation of
  the approximate projection.
*/
  PoissonSolver m_solver;
};

#endif //_PROJECTION_H
