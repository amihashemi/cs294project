#include "DeltaVelocity.H"
#include "FieldData.H"


DeltaVelocity::DeltaVelocity()
{
  // do nothing
}

/* DeltaVelocity::~DeltaVelocity()
{
  // do nothing
  }*/


void DeltaVelocity::init(const FieldData& a_vel){
  m_grid = a_vel.m_grid.grow(-1); 
      // the other stuff doesn't work. fixing it ;-)
      // -Johnny
  m_data.define(a_vel.m_data.getDBox().grow(-1));
  m_data.setVal(0.);
  m_grid = a_vel.m_grid;
}


void DeltaVelocity::operator*=(double a_scalar)
{
  DBox box = m_data.getDBox();
  for(Point pt = box.getLowCorner(); box.notDone(pt); box.increment(pt))
    {
      for (unsigned int comp = 0; comp < DIM; comp++)
      {
        m_data(pt, comp) = m_data(pt, comp) * a_scalar;
      }
    }
}

//const RectMDArray<double, DIM>& DeltaVelocity::operator[](int a_component) const{
//  return m_data[a_component];
//}
//
//RectMDArray<double, DIM>& DeltaVelocity::operator[](int a_component) {
//  return m_data[a_component];
//}


void DeltaVelocity::increment(const double& a_scalar,
			      const DeltaVelocity& a_fieldIncrement){
  DBox bi =a_fieldIncrement.m_data.getDBox();
  RectMDArray<double, DIM> temp(bi);
  a_fieldIncrement.m_data.copyTo(temp);
  for(Point pt = bi.getLowCorner(); bi.notDone(pt); bi.increment(pt)) {
    for (unsigned int comp = 0; comp < DIM; comp++)
    {
      temp(pt, comp) = temp(pt, comp) * a_scalar;
      m_data(pt, comp) = m_data(pt, comp) + temp(pt, comp);
    }
  }
 }

//use the MDArray += and *= in order to do this, it will do all of this on the valid region
