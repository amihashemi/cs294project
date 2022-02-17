// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#include "FieldData.H"
#include "DeltaVelocity.H"

FieldData::FieldData()
{
};

FieldData::FieldData(std::shared_ptr< FFT1D > a_fft1dptr, int a_ghosts)
{
  m_ghosts = a_ghosts;
  m_fft1dPtr = a_fft1dptr;
  m_M = a_fft1dptr->getM();
  m_N = a_fft1dptr->getN();
  int lowerCorner[2] = {0, 0};
  int upperCorner[2] = {m_N - 1, m_N - 1};
  m_grid = DBox(Point(lowerCorner), Point(upperCorner));
  m_data.define(m_grid.grow(a_ghosts));
}


void FieldData::copyTo(FieldData& a_FieldData) const
{
  a_FieldData.m_data.define(m_grid.grow(m_ghosts));
  m_data.copyTo(a_FieldData.m_data);
  a_FieldData.m_M = m_M;
  a_FieldData.m_N = m_N;
  a_FieldData.m_ghosts = m_ghosts;
  a_FieldData.m_grid = DBox(m_grid);
  a_FieldData.m_fft1dPtr = m_fft1dPtr;
}

void FieldData::setBoundaries(RectMDArray< double, DIM >& a_array)
{
  // set Dirichlet BC on all our walls
  Point highCorner = a_array.getDBox().getHighCorner();
  Point lowCorner = a_array.getDBox().getLowCorner();

  for (int dir = 0; dir < DIM; dir++)
  {
    Point unit = getUnitv(dir);
    Point side_a = Point(lowCorner);
    Point side_b = Point(highCorner);

    for (int k = lowCorner[dir]; k <= highCorner[dir]; k++)
    {
      // this loop is to force the BC on all walls for each component
      for (int d = 0; d < DIM; d++)
      {
        a_array(side_a, d) = 0;
        a_array(side_b, d) = 0;
      }
      side_a += unit;
      side_b -= unit;
    }
  }
};

void FieldData::setBoundaries()
{
  setBoundaries(m_data);
};

void FieldData::increment(const double& a_scalar,
                          const DeltaVelocity& a_fieldIncrement)
{
  DBox box = a_fieldIncrement.m_data.getDBox();

  RectMDArray< double, DIM > temp(box);
  for (Point pt = box.getLowCorner(); box.notDone(pt); box.increment(pt))
  {
    for (unsigned int comp = 0; comp < DIM; comp++)
    {

      temp(pt, comp) = a_fieldIncrement.m_data(pt, comp);
      temp(pt, comp) = temp(pt, comp) * a_scalar;
      m_data(pt, comp) = m_data(pt, comp) + temp(pt, comp);
    }
  }
}

