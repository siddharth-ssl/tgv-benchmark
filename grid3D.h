
#ifndef  _GRID_BCC_3D_BASIC_H_
#define  _GRID_BCC_3D_BASIC_H_ 
#include<mm_malloc.h>
#include<iostream>
template <int numfield, typename dataType> 
class gridBCC3D 
 {
  public:
   gridBCC3D(int size1=1, int size2=1, int size3 =1, int stride1=0, int stride2=0, int stride3=0)
  {
   m1 = size1;
   m2 = size2;
   m3 = size3;
   actualGridSize = m1*m2*m3;
   nB1 = stride1;
   nB2 = stride2;
   nB3 = stride3;
   n1 = size1+2*stride1;
   n2 = size2+2*stride2;
   n3 = size3+2*stride3;
   nE1 = size1+ stride1-1;
   nE2 = size2+ stride2-1;
   nE3 = size3+ stride3-1;
   sizeGrid = n1*n2*n3;
   numField = numfield; 
   numElem= sizeGrid*numField ;
   // For vectorization purpose aligned version of memory allocation with vector size as 4096
   data = (dataType*) _mm_malloc(numElem*sizeof(dataType), 4096);
   numField_offset = sizeGrid*numField;
  }
     
  ~gridBCC3D()
  { 
   delete data; 
  }
  
  int size() const
  {
   return numElem;
  }
     
  int getIndex (int i1, int i2, int i3, int dv) const
  { 
   return  ((i3*n2 + i2)*n1 + i1)*numField + dv;
  }
  
  dataType  Node(const int i1,const int i2,const int i3,  int dv) const  { return data[ getIndex(i1,i2,i3,dv)];} 
  dataType& Node(const int i1,const int i2,const int i3,  int dv)        { return data[ getIndex(i1,i2,i3,dv)];} 
   
  dataType value(const int i) const { return data[i];} 
  dataType& value(const int i)      { return data[i];} 
   
  
  
  void initializeIntegers()
  {
   for (int k=nB3; k<=nE3; k++)
   for (int j=nB2; j<=nE2; j++) 
   for (int i=nB1; i<=nE1; i++)
   for (int dv=0; dv<numfield; dv++)
    Node(i,j,k, dv) = getIndex(i,j,k,dv);
  }
   
  template <int N, int bL, typename dataType1>
  friend std::ostream& operator<< (std::ostream& os,const gridBCC3D<N, dataType1> &thisVec);
 
  public:
  int m1, m2, m3;                // Size of physical  grid
  int nB1, nB2, nB3;             // Number of ghost nodes at ends in  direction i and the beginning of the physical grid
  int nE1, nE2, nE3;             // Last element of physical grid nEi = mi + nBi-1
  int n1, n2, n3;                // Size of computational  grid in  direction i, ni = mi+ 2*nBi
  int numElem;                   // total number of element in data array sizeGrid*numField
  int sizeGrid;                  // total number of computational grid points n1*n2*n3
  int numField,actualGridSize; 
  int numField_offset; 
  dataType *data; 
 };
 
	
#endif
