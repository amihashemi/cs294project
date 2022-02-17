#ifndef _ADVECTIONOPERATOR_H_
#define _ADVECTIONOPERATOR_H_
#include <iostream>
#include <cmath>
#include <vector>
#include "DBox.H"
#include "RectMDArray.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
void advectionOperator(DeltaVelocity& a_divuu,
                       const FieldData& a_velocity,
                       const DBox m_grid,
                       const double& a_h);
#endif
