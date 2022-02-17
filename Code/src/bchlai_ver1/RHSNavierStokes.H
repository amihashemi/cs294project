// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#ifndef _RHSNAVIERSTOKES_H
#define _RHSNAVIERSTOKES_H

#include "DeltaVelocity.H"
#include "FieldData.H"

/// Computes the RHS of the Navier-Stokes equation
class RHSNavierStokes
{
public:
  void operator()(DeltaVelocity& a_newDv, const double& a_time,
                  const double& a_dt, const FieldData& a_velocity,
                  DeltaVelocity& a_oldDv);
};

#endif //_RHSNAVIERSTOKES_H
