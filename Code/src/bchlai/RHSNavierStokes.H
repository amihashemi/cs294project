// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#ifndef _RHSNAVIERSTOKES_H
#define _RHSNAVIERSTOKES_H

#include "AdvectionSolver.H"
#include "DeltaVelocity.H"
#include "FieldData.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "Projection.H"
#include "RectMDArray.H"

/// Computes the RHS of the Navier-Stokes equation
class RHSNavierStokes
{
public:
  void operator()(DeltaVelocity& a_k, const double& a_time,
                  const double& a_dt, const FieldData& a_velocity);
};

#endif //_RHSNAVIERSTOKES_H
