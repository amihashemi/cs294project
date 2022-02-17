// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#ifndef _ADVECTIONSOLVER_H_
#define _ADVECTIONSOLVER_H_

#include "DBox.H"
#include "DeltaVelocity.H"
#include "FieldData.H"

// compute the advection term
void advectionSolver(DeltaVelocity& a_divuu, const FieldData& a_velocity,
                     const DBox m_grid, const double& a_h);

#endif
