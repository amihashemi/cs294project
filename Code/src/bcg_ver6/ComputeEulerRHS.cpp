#include "ComputeEulerRHS.H"
#include "FFTW1D.H"
#include "WriteRectMDArray.H"
#include "FieldData.H"

void ComputeEulerRHS::operator() (DeltaVelocity& a_k,
				  const double& a_time, 
				  const double& a_dt,
				  const FieldData& a_velocity) {
  double a_scalar = 1.0;
  
  //std::shared_ptr<FFT1D> fft1dptr(a_velocity.m_N);  std::shared_ptr<FFT1D> fft1dptr(new FFTW1D(a_velocity.m_M));

  //create dummy variables
  FieldData velTemp;
  a_velocity.copyTo(velTemp);
  DeltaVelocity deltaV;
  //int index[DIM];
    double scale = 1./15000.;
  Projection projVel(a_velocity.m_fft1dptr);
  DBox bx = a_velocity.m_grid;
    
  DBox bx0 = bx.grow(-1);
  
  //increment veltemp by a_k
  velTemp.increment(a_scalar,a_k);
    
    velTemp.setBoundaries();
  
  //apply advection on veltemp (what is uu?) advectionOperator(DeltaVelocity& a_divuu,const FieldData& a_velocity,const Box m_grid, const double& a_h);
  advectionOperator(deltaV,velTemp, a_velocity.m_grid,a_dt);
  velTemp.increment(-scale*a_dt,deltaV);
  velTemp.setBoundaries();
  projVel.applyProjection(velTemp);
  velTemp.setBoundaries();
  //assigning a_k :: veltemp2 - a_velocity
  for(Point pt = bx0.getLowCorner(); bx0.notDone(pt); bx0.increment(pt)) {
    /*  for(int i = 0; i <a_velocity.m_M; i++) {
	for(int j = 0; j <a_velocity.m_N; j++) {
	index[0]= i;
	index[1] = j;
	a_k[0][index]= (velTemp2[0])[index] - (a_velocity[0])[index];
	a_k[1][index]= (velTemp2[1])[index] - (a_velocity[1])[index];
	}
	}*/
    a_k.m_data(pt, 0) = (velTemp.m_data(pt, 0) - a_velocity.m_data(pt, 0));
    a_k.m_data(pt, 1) = (velTemp.m_data(pt, 1) - a_velocity.m_data(pt, 1));
  }
  // MDWrite(a_k[1]);
}
