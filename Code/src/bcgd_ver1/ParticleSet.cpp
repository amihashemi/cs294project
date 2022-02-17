// Created by Júlio Caineta on 31/10/17.
// University of Pittsburgh

#include "ParticleSet.H"
#include "../Hockney/Hockney.H"


void ParticleShift::init(const ParticleSet& a_particles)
{
  m_particles.clear();
  m_particles.resize(a_particles.m_particles.size(), DX());
};

void ParticleShift::increment(double a_scale, const ParticleShift& a_rhs)
{
  unsigned int i = 0;
  for (auto p = m_particles.begin(); p != m_particles.end(); ++p, i++)
  {
    (*p).increment(a_scale, a_rhs.m_particles[i]);
  }
};

void ParticleShift::operator*=(double a_scale)
{
  for (auto p = m_particles.begin(); p != m_particles.end(); ++p)
  {
    (*p) *= a_scale;
  }
};

void ParticleShift::setToZero()
{
  std::fill(m_particles.begin(), m_particles.end(), DX());
};

ParticleSet::ParticleSet(shared_ptr< ConvKernel >& a_kerptr, DBox& a_box,
                         double& a_dx, array< double, DIM >& a_lowCorner,
                         int a_M)
{
  m_dx = a_dx;
  m_box = a_box;
  m_lowCorner = a_lowCorner;
  m_hockney.define(a_kerptr, a_dx, a_M);
};

void ParticleSet::increment(const ParticleShift& a_shift)
{
  unsigned int i = 0;
  for (auto p = m_particles.begin(); p != m_particles.end(); ++p, i++)
  {
    (*p).increment(a_shift.m_particles[i]);
  }
};
