// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#include "AdvectionSolver.H"

void advectionSolver(DeltaVelocity& a_divuu,
                     const FieldData& a_velocity,
                     const DBox a_grid,
                     const double& a_h)
{
  int pt[DIM];
  int downShift[DIM];
  int upShift[DIM];
  FieldData velocity;
  a_velocity.copyTo(velocity);

  a_divuu.init(a_velocity);
  for (int i = 0; i < a_grid.size(0); i++)
  {
    pt[0] = i;
    for (int j = 0; j < a_grid.size(1); j++)
    {
      pt[1] = j;
      for (int k = 0; k < DIM; k++)
      {
        for (int l = 0; l < DIM; l++)
        {
          double left;
          double right;
          double diff;

          upShift[0] = pt[0];
          upShift[1] = pt[1];
          downShift[0] = pt[0];
          downShift[1] = pt[1];
          upShift[l] += 1;
          downShift[l] -= 1;

          left = (velocity.m_data(Point(pt), k)
                  + velocity.m_data(Point(upShift), k))
                 * (a_velocity.m_data(Point(pt), l)
                    + a_velocity.m_data(Point(upShift), l)) / 4;
          right = (velocity.m_data(Point(pt), k)
                   + velocity.m_data(Point(downShift), k))
                  * (velocity.m_data(Point(pt), l)
                     + velocity.m_data(Point(downShift), l)) / 4;

          diff = (left - right) / a_h;

          a_divuu.m_data(Point(pt), k) =
              a_divuu.m_data(Point(pt), k) + diff;
        }
      }
    }
  }
}
