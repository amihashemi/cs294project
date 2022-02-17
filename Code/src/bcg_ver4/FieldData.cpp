#include "FieldData.H"
#include "DeltaVelocity.H"

FieldData::FieldData(){
}

//FIELD DATA is Q from the notes
FieldData::FieldData(std::shared_ptr<FFT1D> a_fft1dptr, int a_nghosts){
  m_ghosts = a_nghosts;
  m_fft1dptr = a_fft1dptr;
  m_M = a_fft1dptr->getM();
  m_N = a_fft1dptr->getN();
  int lowerCorner[2] = {0,0};
  int upperCorner[2] = {m_N-1, m_N-1};
  m_grid = DBox(Point(lowerCorner), Point(upperCorner));
      m_data.define(m_grid.grow(a_nghosts));
//make a box with low corner = 00, high corner = m_n, m_M?? make this box m_grid
  //  m_grid = m_data.getDBox();
}

//const RectMDArray<double, DIM>& FieldData::operator[](const Point a_component) const{
//  return m_data[a_component];
//}

//RectMDArray<double, DIM>& FieldData::operator[](const Point a_component){
//  return m_data[a_component];
//}

void FieldData::copyTo(FieldData& a_FieldData) const
{
      a_FieldData.m_data.define(m_grid.grow(m_ghosts));
      m_data.copyTo(a_FieldData.m_data);
  a_FieldData.m_M = m_M;
  a_FieldData.m_N = m_N;
  a_FieldData.m_ghosts = m_ghosts;
  a_FieldData.m_grid = DBox(m_grid);
  a_FieldData.m_fft1dptr = m_fft1dptr;
}

//void FieldData::fillGhosts(RectMDArray<double, DIM>& array){
////  int highCorner[DIM];
////  int lowCorner[DIM];
//  Point highCorner = array.getDBox().getHighCorner();
//  Point lowCorner = array.getDBox().getLowCorner();
//  int ghostx = lowCorner[0];
//  int ghosty = lowCorner[1];
//  //fill Ghosts for bottom row, top row
//  for(ghostx=lowCorner[0]+1; ghostx <highCorner[0]; ghostx++){
//    int bottomIndex[DIM] = {ghostx, highCorner[1]-1};
//    int bottomGhostIndex[DIM] = {ghostx, lowCorner[1]};
//    int topIndex[DIM] = {ghostx, lowCorner[1] + 1};
//    int topGhostIndex[DIM] = {ghostx, highCorner[1]};
//    array[Point(bottomGhostIndex)];
//    array[Point(bottomGhostIndex)] = array[Point(bottomIndex)];
//    array[Point(topGhostIndex)] = array[Point(topIndex)];
//  }
//  // fillGhosts for left row, right row
//  for(ghosty=lowCorner[1]+1; ghosty < highCorner[1]; ghosty++){
//    int leftIndex[DIM] = {highCorner[0] -1, ghosty};
//    int leftGhostIndex[DIM] = {lowCorner[0], ghosty};
//    int rightIndex[DIM] = {lowCorner[0]+1, ghosty};
//    int rightGhostIndex[DIM] = {highCorner[0], ghosty};
//    array[Point(leftGhostIndex)] = array[Point(leftIndex)];
//    array[Point(rightGhostIndex)] = array[Point(rightIndex)];
//    }
//}
//
//void FieldData::fillGhosts(){
//  fillGhosts(m_data);
//  // DeltaVelocities do not have the ghost data
//  // Fill in the edges from the mirror on the opposite side of the box
//  //Does this mean that I want to expand m_grid to have extra spaces around the outside?
//}


void FieldData::increment(const double& a_scalar,
			  const DeltaVelocity& a_fieldIncrement){
  DBox bi =a_fieldIncrement.m_data.getDBox();
  RectMDArray<double, DIM> temp(bi);
  for(Point pt = bi.getLowCorner(); bi.notDone(pt); bi.increment(pt)) {
    for (unsigned int comp = 0; comp < DIM; comp++)
    {
      temp(pt, comp) = temp(pt, comp) * a_scalar;
      m_data(pt, comp) = m_data(pt, comp) + temp(pt, comp);
    }
  }
}

