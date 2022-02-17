// Final project for CS 294-73 at Berkeley
// Amirreza Hashemi and Júlio Caineta
// 2017

#include "FieldData.H"


FieldData::FieldData(){
}


FieldData::FieldData(DBox a_grid, int a_nComponent, int a_ghost, int a_M, int a_N){
    m_ghosts = a_nghosts;
    m_M = a_M;
    m_N = a_N;
    int lowerCorner[2] = {0,0};
    int upperCorner[2] = {m_N-1, m_N-1};
    m_grid = Box(lowerCorner, upperCorner);
    for (int i = 0; i<DIM; i++)
    {
        m_data[i].define(m_grid.grow(a_nghosts));
    }

}

void FieldData::increment(const double& a_scalar,
                              const DeltaVelocity& a_fieldIncrement){
    
    DBox bi =a_fieldIncrement.m_data[0].getBox();
    RectMDArray <double, bi> temp;
    for (int i = 0; i <DIM; i++){
        a_fieldIncrement.m_data[i].copyTo(temp);
        temp*=a_scalar;
        m_data[i] += temp;
        
    }
}


