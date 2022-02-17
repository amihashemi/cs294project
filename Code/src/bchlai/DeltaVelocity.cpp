// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#include "DeltaVelocity.H"


DeltaVelocity::DeltaVelocity()
{};

void DeltaVelocity::init(const FieldData& a_velocity)
{
  m_grid = a_velocity.m_grid.grow(-1);
  m_data.define(a_velocity.m_data.getDBox().grow(-1));
  m_data.setVal(0.);
  m_grid = a_velocity.m_grid;
};


void DeltaVelocity::operator*=(double a_scalar)
{
  DBox box = m_data.getDBox();
  for (Point pt = box.getLowCorner(); box.notDone(pt); box.increment(pt))
  {
    for (unsigned int comp = 0; comp < DIM; comp++)
    {
      m_data(pt, comp) = m_data(pt, comp) * a_scalar;
    }
  }
};


void DeltaVelocity::increment(const double& a_scalar,
                              const DeltaVelocity& a_fieldIncrement)
{
  DBox box = a_fieldIncrement.m_data.getDBox();

  for (Point pt = box.getLowCorner(); box.notDone(pt); box.increment(pt))
  {
    for (unsigned int comp = 0; comp < DIM; comp++)
    {
      m_data(pt, comp) = m_data(pt, comp)
                         + a_fieldIncrement.m_data(pt, comp) * a_scalar;
    }
  }
};
