// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#include "DeltaVelocity.H"

DeltaVelocity::DeltaVelocity()
{
    
}


void DeltaVelocity::increment(const double& a_scalar,
                              const DeltaVelocity& a_fieldIncrement){
    
    DBox bi =a_fieldIncrement.m_data[0].getBox();
    RectMDArray <double, bi> temp;
    for (int i = 0; i <DIM; i++){
        a_fieldIncrement.m_data[i].copyTo(temp);
        temp*=a_scalar;
        m_data[i] += temp;
        
    }
}
