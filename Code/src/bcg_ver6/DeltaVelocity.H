#ifndef _DELTAVELOCITY_H_
#define _DELTAVELOCITY_H_

#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "DBox.H"
#include "RectMDArray.H"
class FieldData;
using namespace std;
class DeltaVelocity
{
public:
  DeltaVelocity();
  //~DeltaVelocity();
  void init(const FieldData& a_vel);

  // This will give you X or Y component
//  const RectMDArray<double, DIM>& operator[](int a_component) const;
//  RectMDArray<double, DIM>& operator[](int a_component);
  void increment(const double& a_scalar,
                 const DeltaVelocity& a_fieldIncrement);
  void operator*=(double a_scalar);

  DBox m_grid;

  RectMDArray<double, DIM> m_data;

private:
  DeltaVelocity(const DeltaVelocity& a_DeltaV);
  void operator=(const DeltaVelocity& a_DeltaV);


};
#endif
