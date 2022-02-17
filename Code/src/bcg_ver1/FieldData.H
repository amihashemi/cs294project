#ifndef _FIELDDATA_H_
#define _FIELDDATA_H_

#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <list>
#include <memory>
#include "DIM.H"
#include "DBox.H"
#include "RectMDArray.H"
#include "FFT1D.H"
using namespace std;
class DeltaVelocity;
class FieldData
{
public:
  FieldData();
  FieldData(std::shared_ptr<FFT1D> a_fft1dptr, int a_nghosts);
  //~FieldData();
//  const RectMDArray<double, DIM>& operator[](const Point a_component) const;
//  RectMDArray<double, DIM>& operator[](const Point a_component);
//  void fillGhosts(RectMDArray<double>& array);
//  void fillGhosts();
  void increment(const double& a_scalar,
                 const DeltaVelocity& a_fieldIncrement);
  void copyTo(FieldData& a_FieldData) const;
  DBox m_grid;
  int m_M,m_N,m_ghosts;
  RectMDArray<double, DIM> m_data;
  std::shared_ptr<FFT1D> m_fft1dptr;
  
private:
  FieldData(const FieldData& a_FieldData);
  void operator=(const FieldData& a_FieldData);

};
#endif
