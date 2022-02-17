// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#ifndef _DELTAVELOCITY_H_
#define _DELTAVELOCITY_H_

#include "FieldData.H"
#include "RectMDArray.H"

class FieldData;

class DeltaVelocity
{
public:
  //Default constructor
  DeltaVelocity();

  //constructor
  DeltaVelocity(DBox a_grid);

  //destructor
  ~DeltaVelocity()
  {};

  RectMDArray< double, DIM >& getVelocity();

  void init(const FieldData& a_velocity);

  void increment(const double& a_scalar,
                 const DeltaVelocity& a_fieldIncrement);

  void operator*=(double a_scalar);

  DBox m_grid;
  RectMDArray< double, DIM > m_data;

};

#endif
