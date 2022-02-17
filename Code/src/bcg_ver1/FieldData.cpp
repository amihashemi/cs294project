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
  m_grid = DBox(lowerCorner, upperCorner);
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
//
//void FieldData::fillGhosts(RectMDArray<double>& array){
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
//    RectMDArray<double> temp;
//    temp =  m_data[1];
//  fillGhosts(temp);
//   // fillGhosts(m_data[1]);
//  // DeltaVelocities do not have the ghost data
//  // Fill in the edges from the mirror on the opposite side of the box
//  //Does this mean that I want to expand m_grid to have extra spaces around the outside?
//}
////
//
//
//void FieldData::fillGhosts(RectMDArray<double>& array){
//    //  int highCorner[DIM];
//    //  int lowCorner[DIM];
//    Point highCorner = array.getDBox().getHighCorner();
//    Point lowCorner = array.getDBox().getLowCorner();
//    int ghostx = lowCorner[0];
//    int ghosty = lowCorner[1];
//    //fill Ghosts for bottom row, top row
//    for(ghostx=lowCorner[0]+1; ghostx <highCorner[0]; ghostx++){
//        int bottomIndex[DIM] = {ghostx, highCorner[1]-1};
//        int bottomGhostIndex[DIM] = {ghostx, lowCorner[1]};
//        int topIndex[DIM] = {ghostx, lowCorner[1] + 1};
//        int topGhostIndex[DIM] = {ghostx, highCorner[1]};
//        array[Point(bottomGhostIndex)] = array[Point(bottomIndex)];
//        array[Point(topGhostIndex)] = array[Point(topIndex)];
//    }
//    // fillGhosts for left row, right row
//    for(ghosty=lowCorner[1]+1; ghosty < highCorner[1]; ghosty++){
//        int leftIndex[DIM] = {highCorner[0] -1, ghosty};
//        int leftGhostIndex[DIM] = {lowCorner[0], ghosty};
//        int rightIndex[DIM] = {lowCorner[0]+1, ghosty};
//        int rightGhostIndex[DIM] = {highCorner[0], ghosty};
//        array[Point(leftGhostIndex)] = array[Point(leftIndex)];
//        array[Point(rightGhostIndex)] = array[Point(rightIndex)];
//    }
//}

//void FieldData::fillGhosts(){
////    RectMDArray<double> temp;
////    temp =  m_data[1];
////    fillGhosts(temp);
//    Point highCorner = array.getDBox().getHighCorner();
//    Point lowCorner = array.getDBox().getLowCorner();
//    int ghostx = lowCorner[0];
//    int ghosty = lowCorner[1];
//    //fill Ghosts for bottom row, top row
//    for(ghostx=lowCorner[0]+1; ghostx <highCorner[0]; ghostx++){
//        int bottomIndex[DIM] = {ghostx, highCorner[1]-1};
//        int bottomGhostIndex[DIM] = {ghostx, lowCorner[1]};
//        int topIndex[DIM] = {ghostx, lowCorner[1] + 1};
//        int topGhostIndex[DIM] = {ghostx, highCorner[1]};
//        m_data[Point(bottomGhostIndex)][0] = m_data[Point(bottomIndex)][0];
//        m_data[Point(topGhostIndex)][0] = m_data[Point(topIndex)][0];
//    }
//    // fillGhosts for left row, right row
//    for(ghosty=lowCorner[1]+1; ghosty < highCorner[1]; ghosty++){
//        int leftIndex[DIM] = {highCorner[0] -1, ghosty};
//        int leftGhostIndex[DIM] = {lowCorner[0], ghosty};
//        int rightIndex[DIM] = {lowCorner[0]+1, ghosty};
//        int rightGhostIndex[DIM] = {highCorner[0], ghosty};
//        m_data[0][Point(leftGhostIndex)] = m_data[0][Point(leftIndex)];
//        m_data[0][Point(rightGhostIndex)] = m_data[0][Point(rightIndex)];
//    }
//}


void FieldData::increment(const double& a_scalar,
			  const DeltaVelocity& a_fieldIncrement){
  DBox bi =a_fieldIncrement.m_data.getDBox();
  RectMDArray<double> temp(bi);
  for(Point pt = bi.getLowCorner(); bi.notDone(pt); bi.increment(pt)) {
    temp[pt] = temp[pt] * a_scalar;
    m_data[pt] = m_data[pt] + temp[pt];
  }
}

// (ParticleShift& a_kOut, 
//                      const double& a_time, const double& dt,
//                      const ParticleSet& a_state,
//                      const ParticleShift& a_kIn)

