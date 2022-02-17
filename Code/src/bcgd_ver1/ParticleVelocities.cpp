// Created by Júlio Caineta on 31/10/17.
// University of Pittsburgh
#include "ParticleVelocities.H"

ParticleVelocities::ParticleVelocities()
{};

void ParticleVelocities::operator()(ParticleShift& a_k, const double& a_time,
                                    const double& dt, ParticleSet& a_state)
{
  double h = a_state.m_dx;
  array< int, 2 > i_x = {{1, 0}};
  array< int, 2 > i_y = {{0, 1}};
  array< int, 2 > i_xy = {{1, 1}};

  Point pt_x = Point(i_x);     // (1, 0)
  Point pt_y = Point(i_y);     // (0, 1)
  Point pt_xy = Point(i_xy);   // (1, 1)

  // 1. depositing the charges in the particles on the grid
  vector< Particle > particles(a_state.m_particles);
  RectMDArray< double > omega(a_state.m_box);
  omega.setVal(0.);

  unsigned int k = 0;

  for (auto p = particles.begin(); p != particles.end(); ++p, k++)
  {
    (*p).increment(a_k.m_particles[k]);  // x_k
    array< int, 2 > i_k;
    array< int, 2 > s_k;

    for (unsigned int d = 0; d < 2; d++)
    {
      i_k[d] = std::floor(((*p).m_x[d] - a_state.m_lowCorner[d]) / h);
      s_k[d] = ((*p).m_x[d] - i_k[d] * h) / h;
    }

    Point pt = Point(i_k);

    omega[pt] += (*p).strength * (1 - s_k[0]) * (1 - s_k[1]);
    omega[pt + pt_x] += (*p).strength * s_k[0] * (1 - s_k[1]);
    omega[pt + pt_y] += (*p).strength * (1 - s_k[0]) * s_k[1];
    omega[pt + pt_xy] += (*p).strength * s_k[0] * s_k[1];
  }

  // 2. convolution with Green's function, using Hockney's method
  a_state.m_hockney.convolve(omega);

  // 3. compute the velocity fields on the grid using finite differences
  DBox D(a_state.m_box.getLowCorner(), a_state.m_box.getHighCorner());
//  RectMDArray< double, 2 > ug(D);
  RectMDArray< double> ugx(D);
  RectMDArray< double> ugy(D);
//  ug.setVal(0.);
  ugx.setVal(0.);
  ugy.setVal(0.);
  DBox domain(D);
  domain = domain.grow(-1);

  for (Point pt = domain.getLowCorner();
       domain.notDone(pt); domain.increment(pt))
  {
//    ug(pt, 1) = (omega[pt + pt_y] - omega[pt - pt_y]) / (2 * h);
//    ug(pt, 2) = (omega[pt - pt_x] - omega[pt + pt_x]) / (2 * h);
    ugx[pt] = (omega[pt + pt_y] - omega[pt - pt_y]) / (2 * h);
    ugy[pt] = (omega[pt - pt_x] - omega[pt + pt_x]) / (2 * h);
  }

  // 4. interpolate the fields from the grids to the particles
  k = 0;
  for (auto p = particles.begin(); p != particles.end(); ++p, k++)
  {
    array< int, 2 > i_k;
    array< int, 2 > s_k;

    for (unsigned int d = 0; d < 2; d++)
    {
      i_k[d] = std::floor(((*p).m_x[d] - a_state.m_lowCorner[d]) / h);
      s_k[d] = ((*p).m_x[d] - i_k[d] * h) / h;
    }

    Point pt = Point(i_k);

//    for (unsigned int d = 1; d < 3; d++)
//    {
//      double velo = ug(pt, d) * (1 - s_k[0]) * (1 - s_k[1])
//                    + ug(pt + pt_x, d) * s_k[0] * (1 - s_k[1])
//                    + ug(pt + pt_y, d) * (1 - s_k[0]) * s_k[1]
//                    + ug(pt + pt_xy, d) * s_k[0] * s_k[1];
//      a_k.m_particles[k].m_x[d - 1] = dt * velo;
//
//      if (std::abs(dt * velo) > 0.0000001)
//      {
//        std::cout << "outch!!" << endl;
//        std::cin.get();
//      }
//    }

    double velox = ugx[pt] * (1 - s_k[0]) * (1 - s_k[1])
                  + ugx[pt + pt_x] * s_k[0] * (1 - s_k[1])
                  + ugx[pt + pt_y] * (1 - s_k[0]) * s_k[1]
                  + ugx[pt + pt_xy] * s_k[0] * s_k[1];

    a_k.m_particles[k].m_x[0] = dt * velox;
    
    double veloy = ugy[pt] * (1 - s_k[0]) * (1 - s_k[1])
                  + ugy[pt + pt_x] * s_k[0] * (1 - s_k[1])
                  + ugy[pt + pt_y] * (1 - s_k[0]) * s_k[1]
                  + ugy[pt + pt_xy] * s_k[0] * s_k[1];

    a_k.m_particles[k].m_x[1] = dt * veloy;

  }
};