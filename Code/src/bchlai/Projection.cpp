// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#include "RectMDArray.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "FieldData.H"
#include "DeltaVelocity.H"
#include "PoissonSolver.H"
#include "Projection.H"
#include "WriteRectMDArray.H"

using namespace std;

Projection::Projection(std::shared_ptr< FFT1D > a_fft1dPtr)
{
  m_solver = PoissonSolver(a_fft1dPtr);
  m_N = a_fft1dPtr->getN();
}

void Projection::applyProjection(FieldData& a_velocity) const
{
  DBox vBox(a_velocity.m_grid);
  RectMDArray< double > R(vBox);
  divergence(R, a_velocity);

  m_solver.solve(R);

  RectMDArray< double > Rghost(vBox.grow(1));
  R.copyTo(Rghost);
  DeltaVelocity dv;
  dv.init(a_velocity);
  gradient(dv, Rghost);
  a_velocity.increment(-1.0, dv);
}

void Projection::gradient(DeltaVelocity& a_vector,
                          const RectMDArray< double >& a_scalar) const
{
  Point lowCorner = a_vector.m_grid.getLowCorner();
  Point highCorner = a_vector.m_grid.getHighCorner();
  double h = 1.0 / m_N;
  DBox bx = a_vector.m_grid;

  for (Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt))
  {
    int i = pt[0];
    int j = pt[1];
    int plusE0[2] = {i + 1, j};  // i+e^0 (right side)
    int minusE0[2] = {i - 1, j}; // i-e^0 (left side)
    int plusE1[2] = {i, j + 1};  // i+e^1 (top side)
    int minusE1[2] = {i, j - 1}; // i-e^1 (bottom side)

    a_vector.m_data(pt, 0) =
        (a_scalar[Point(plusE0)] - a_scalar[Point(minusE0)]) / (2 * h);
    a_vector.m_data(pt, 1) =
        (a_scalar[Point(plusE1)] - a_scalar[Point(minusE1)]) / (2 * h);
  }
}

void Projection::divergence(RectMDArray< double >& a_scalar,
                            const FieldData& a_vector) const
{
  DBox bx = a_vector.m_grid;
  double h = 1.0 / a_vector.m_N;

  for (Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt))
  {
    int i = pt[0];
    int j = pt[1];
    int plusE0[2] = {i + 1, j};  // i+e^0 (right side)
    int minusE0[2] = {i - 1, j}; // i-e^0 (left side)
    int plusE1[2] = {i, j + 1};  // i+e^1 (top side)
    int minusE1[2] = {i, j - 1}; // i-e^1 (bottom side)

    double xPortion =
        a_vector.m_data(Point(plusE0), 0) - a_vector.m_data(Point(minusE0), 0);
    double yPortion =
        a_vector.m_data(Point(plusE1), 1) - a_vector.m_data(Point(minusE1), 1);
    a_scalar[pt] = (xPortion + yPortion) / (2 * h);
  }
}
