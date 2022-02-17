#include <cstdio>
#include <iostream>
#include "WriteRectMDArray.H"
#include "DBox.H"
#include "RectMDArray.H"
#include "VisitWriter.H"
using namespace std;


const char* MDWrite(RectMDArray<double,1>* a_array)
{
  static int fileCount = 0;
  static char nameBuffer[10];

  if (a_array == NULL)
    {
      return nameBuffer;
    }

  sprintf(nameBuffer, "md%06d",fileCount);
  MDWrite(nameBuffer, a_array);

  fileCount++;

  return nameBuffer;
};

void MDWrite(const char*           a_filename,
             RectMDArray<double>* a_array)
{
  FILE* fp = vtk_open_file(a_filename);
  double origin[3]={0,0,0};
  double dx=1.0;
  char* vars[1];
  vars[0] = "0";
  MDWrite(fp, *a_array, vars, origin, dx);
  }
