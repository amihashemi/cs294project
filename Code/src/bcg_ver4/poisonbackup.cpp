  
#include "PoissonSolver.H"

PoissonSolver::PoissonSolver()
{
}
PoissonSolver::PoissonSolver(std::shared_ptr<FFT1D> a_fft1dPtr)
{
  
  m_fft1dptr = a_fft1dPtr;
  m_N = a_fft1dPtr->getN();
  m_M = a_fft1dPtr->getM();
}

void PoissonSolver::solve(RectMDArray<double>& a_Rhs) const
{
	//create a contructor
  DBox bx = a_Rhs.getDBox();
  RectMDArray<complex<double> > fftwForward(bx);
  RectMDArray<complex<double> > fourierCoef(bx);
  
  double h = 1/m_M;

	for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt))
//  for (int k = 0; k < bx.sizeOf();k++)
  {
		fftwForward[pt] = complex<double>(a_Rhs[pt],0.);
  }
  FFTMD fftmd = FFTMD(m_fft1dptr);
  fftmd.forwardCC(fftwForward);
 
	 for (int i= 0; i<m_M; i++)
	 {
		 for (int j= 0; j<m_M; j++)
		 {
			 int index[2] ={i,j};
			 
			 if(i ==0 && j ==0)
			 {
				 fourierCoef[index] =complex<double>(0,0);
			 }
			 else
			 {
				 complex <double> div(2*cos(2*M_PI*i*h) - 2*cos(2*M_PI*j*h - 4),0);
				 fourierCoef[index] = fftwForward[index]/div;
			 }

		 }
	 }


	 fftmd.inverseCC(fourierCoef);

	for(Point pt = bx.getLowCorner(); bx.notDone(pt); bx.increment(pt))
//	 for (int k= 0; k < bx.sizeOf(); k++)
	 {
		
           a_Rhs[pt] = real(fourierCoef[pt])/m_N;
		 
	 }

}

///  int m_M,m_N; DBox m_grid;std::shared_ptr<FFT1D> m_fft1dptr;
