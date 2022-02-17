
#include "WriteMDArray.H"
#include "MDArray.H"
#include "VisitWriter.H"
//#include "ParticleSet.H"

#include <cstdio>

const char* MDWrite(MDArray<float>* a_array)
{
  static int fileCount = 0;
  static char nameBuffer[10];
  if(a_array == NULL)
    {
      return nameBuffer;
    }
  sprintf(nameBuffer, "md%d",fileCount);
  MDWrite(nameBuffer, a_array);
  fileCount++;
  return nameBuffer;
}

/*
const char* PWrite(const ParticleSet* a_array)
{
  static int fileCount = 0;
  static char nameBuffer[10];
  if(a_array == NULL)
    {
      return nameBuffer;
    }
  sprintf(nameBuffer, "PART.%d",fileCount);
  PWrite(nameBuffer, a_array);
  fileCount++;
  return nameBuffer;
}

void PWrite(const char* a_filename, const ParticleSet* a_p)
{
  if(a_filename == NULL || a_p == NULL)
    {
      return;
    }
  unsigned int size = a_p->m_particles.size();
  std::vector<float> x(3*size);
  for(unsigned int i=0; i<size; i++)
    {
      const Particle& p = a_p->m_particles[i];
      x[i*3] = p.m_x[0];
      x[i*3+1] = p.m_x[1];
#if DIM==3
      x[i*3+2] = p.m_x[2];
#else
      x[i*3+2] = 0.0;
#endif
    }

  write_point_mesh(a_filename, 0, size,
		   &(x[0]), 0, 0,
		   0, 0);
}
*/

void MDWrite(const char* a_filename, MDArray<float>* a_array)
{
  if(a_filename == NULL || a_array == NULL)
    {
      return;
    }
  int dim[3] = {1,1,1};
  int vardims[1] ={1};
  int centering[1]={1};
  double* vars[1];
 
  const char * const varnames[] = { "cellCentered" };
  int lo[DIM];
  int hi[DIM];
  const DBox& box = a_array->getBox();
  box.getLowCorner(lo);
  box.getHighCorner(hi);
  for(int i=0; i<DIM;i++)
    {
      dim[i] = hi[i]-lo[i]+1;
    }
  float& val = a_array->operator[](0);
  vars[0] = reinterpret_cast<double*>(&val);
  write_regular_mesh(a_filename, 1, dim, 1, vardims, centering,  varnames, vars);

}
void MDWrite(const char* filename, MDArray<double>& a_array)
{
  DBox bx = a_array.getBox();
  MDArray<float> array(bx);
  for (int k = 0; k < bx.sizeOf();k++)
    {
      array[k] = (float) a_array[k];
    }
  MDWrite(filename,&array);
}
const char* MDWrite(MDArray<double>& a_array)
{
  DBox bx = a_array.getBox();
  MDArray<float> array(bx);
  for (int k = 0; k < bx.sizeOf();k++)
    {
      array[k] = (float) a_array[k];
    }
  return MDWrite(&array);
  
}

