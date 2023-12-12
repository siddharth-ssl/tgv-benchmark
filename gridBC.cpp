template <int N,typename myReal>
void periodicAll(gridBCC3D<N,myReal> &myGrid)
{
 int i2= myGrid.nB2-1;
 for (int i1=0; i1< myGrid.n1;  i1++)
 {
  for (  int i3=0; i3<  myGrid.n3;  i3++)
  {
   for(int k=0;k<N;k++)
   {
    myGrid.Node(i1,i2,i3,k)            = myGrid.Node(i1,myGrid.nE2,i3,k);
    myGrid.Node(i1,myGrid.nE2+1,i3,k)  = myGrid.Node(i1,myGrid.nB2,i3,k);
   }
  }
 }
 int i1 =myGrid.nB1-1;
 for (int i3=0; i3<  myGrid.n3;  i3++)
 {
  for (  int i2=0; i2<  myGrid.n2;  i2++)
  {
   for(int k=0;k<N;k++)
   {
    myGrid.Node(i1,i2,i3,k)            = myGrid.Node(myGrid.nE1,i2,i3,k);
    myGrid.Node(myGrid.nE1+1,i2,i3,k)  = myGrid.Node(myGrid.nB1,i2,i3,k);
   }
  }
 }
 int i3= myGrid.nB3-1;
 for (int i1=0; i1< myGrid.n1;  i1++)
 {
  for (  int i2=0; i2<  myGrid.n2;  i2++)
  {
   for(int k=0;k<N;k++)
   {
    myGrid.Node(i1,i2,i3,k)            = myGrid.Node(i1,i2,myGrid.nE3,k);
    myGrid.Node(i1,i2,myGrid.nE3+1,k)  = myGrid.Node(i1,i2,myGrid.nB3,k);
   }
  }
 }
}