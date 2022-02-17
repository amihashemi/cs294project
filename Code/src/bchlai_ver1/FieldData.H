// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#ifndef _FIELDDATA_H
#define _FIELDDATA_H

#include "DeltaVelocity.H"

class FieldData
{
public:
  // Default constructor
  FieldData();

  // Constructor
  FieldData(DBox a_grid, int a_nComponent, int a_ghost, int a_M, int a_N);

  // Destructor
  ~FieldData();

  void fillGhosts();

  void increment(const double& a_scalar, const DeltaVelocity& a_fieldIncrement);

  void imposeConstraint();

  int m_components;
  DBox m_grid;
  int m_M;
  int m_N;
  int m_ghosts;
  RectMDArray< double, DIM > m_data;
};

#endif //_FIELDDATA_H
