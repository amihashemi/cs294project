// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#include "RHSNavierStokes.H"

void RHSNavierStokes::operator()(DeltaVelocity& a_k,
                                 const double& a_time,
                                 const double& a_dt,
                                 const FieldData& a_velocity)
{
  double a_scalar = 1.0;

  FieldData velTemp;
  a_velocity.copyTo(velTemp);
  DeltaVelocity deltaV;

  double scale = 1. / 40000.;
  Projection projVel(a_velocity.m_fft1dPtr);
  DBox bx = a_velocity.m_grid;

  DBox bx0 = bx.grow(-1);

  velTemp.increment(a_scalar, a_k);

  velTemp.setBoundaries();

  advectionSolver(deltaV, velTemp, a_velocity.m_grid, a_dt);
  velTemp.increment(-scale * a_dt, deltaV);
  velTemp.setBoundaries();
  projVel.applyProjection(velTemp);
  velTemp.setBoundaries();
  for (Point pt = bx0.getLowCorner(); bx0.notDone(pt); bx0.increment(pt))
  {
    a_k.m_data(pt, 0) = velTemp.m_data(pt, 0) - a_velocity.m_data(pt, 0);
    a_k.m_data(pt, 1) = velTemp.m_data(pt, 1) - a_velocity.m_data(pt, 1);
  }
}
