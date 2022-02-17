// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#ifndef _FIELDDATA_H
#define _FIELDDATA_H

#include "DeltaVelocity.H"
#include "FFT1D.H"
#include "RectMDArray.H"

class DeltaVelocity;

class FieldData
{
public:
  // Default constructor
  FieldData();

  // Constructor
  FieldData(DBox a_grid, int a_nComponent, int a_ghost, int a_M, int a_N);

  FieldData(std::shared_ptr< FFT1D > a_fft1dptr, int a_nghosts);

  // Destructor
  ~FieldData()
  {};

  void setBoundaries(RectMDArray< double, DIM >& array);

  void setBoundaries();

  void increment(const double& a_scalar, const DeltaVelocity& a_fieldIncrement);

  void copyTo(FieldData& a_FieldData) const;

  int m_components;
  DBox m_grid;
  int m_M;
  int m_N;
  int m_ghosts;
  RectMDArray< double, DIM > m_data;
  std::shared_ptr< FFT1D > m_fft1dPtr;
};

#endif //_FIELDDATA_H
