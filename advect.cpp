template <int numfield, typename dataType> 
void advection(gridBCC3D< numfield,dataType>  &Grid1,gridBCC3D< numfield,dataType>  &Grid2)
{
 int i1, i2, i3; 
 for(i3 = Grid1.nB3; i3<= Grid1.nE3;  i3++)
  for(i2 = Grid1.nB2; i2<= Grid1.nE2;  i2++)
   for(i1 = Grid1.nB1; i1<= Grid1.nE1;  i1++)
    for(int dv=0;dv<27;dv++)
     Grid2.Node(i1,i2,i3,dv) = Grid1.Node(i1-(int)(cx[dv]),i2-(int)(cy[dv]),i3-((int)cz[dv]),dv); 
}

template <int numfield, typename dataType> 
void advection(gridBCC3D< numfield,dataType>  &Grid1)
{
 int i1, i2, i3; 
for(i3 = Grid1.nE3; i3>= Grid1.nB3;  i3--)
    for(i2 = Grid1.nE2; i2>= Grid1.nB2;  i2--)
      for(i1 = Grid1.nE1; i1>= Grid1.nB1;  i1--) {
        Grid1.Node(i1, i2, i3, DV_P1_ZERO_ZERO) = Grid1.Node(i1-1, i2, i3, DV_P1_ZERO_ZERO);
        Grid1.Node(i1, i2, i3, DV_ZERO_P1_ZERO) = Grid1.Node(i1, i2-1, i3, DV_ZERO_P1_ZERO);
        Grid1.Node(i1, i2, i3, DV_ZERO_ZERO_P1) = Grid1.Node(i1, i2, i3-1, DV_ZERO_ZERO_P1);


        Grid1.Node(i1, i2, i3, DV_P1_P1_ZERO) = Grid1.Node(i1-1, i2-1, i3, DV_P1_P1_ZERO);
        Grid1.Node(i1, i2, i3, DV_M1_P1_ZERO) = Grid1.Node(i1+1, i2-1, i3, DV_M1_P1_ZERO);
        Grid1.Node(i1, i2, i3, DV_ZERO_P1_P1) = Grid1.Node(i1, i2-1, i3-1, DV_ZERO_P1_P1);
        Grid1.Node(i1, i2, i3, DV_ZERO_M1_P1) = Grid1.Node(i1, i2+1, i3-1, DV_ZERO_M1_P1);
        Grid1.Node(i1, i2, i3, DV_P1_ZERO_P1) = Grid1.Node(i1-1, i2, i3-1, DV_P1_ZERO_P1);
        Grid1.Node(i1, i2, i3, DV_M1_ZERO_P1) = Grid1.Node(i1+1, i2, i3-1, DV_M1_ZERO_P1);

        Grid1.Node(i1, i2, i3, DV_P1_P1_P1) = Grid1.Node(i1-1, i2-1, i3-1, DV_P1_P1_P1);
        Grid1.Node(i1, i2, i3, DV_M1_P1_P1) = Grid1.Node(i1+1, i2-1, i3-1, DV_M1_P1_P1);
        Grid1.Node(i1, i2, i3, DV_M1_M1_P1) = Grid1.Node(i1+1, i2+1, i3-1, DV_M1_M1_P1);
        Grid1.Node(i1, i2, i3, DV_P1_M1_P1) = Grid1.Node(i1-1, i2+1, i3-1, DV_P1_M1_P1);
      }

  for(i3 = Grid1.nB3; i3<= Grid1.nE3;  i3++)
    for(i2 = Grid1.nB2; i2<= Grid1.nE2;  i2++)
      for(i1 = Grid1.nB1; i1<= Grid1.nE1;  i1++) {
        Grid1.Node(i1, i2, i3, DV_M1_ZERO_ZERO) = Grid1.Node(i1+1, i2, i3, DV_M1_ZERO_ZERO);
        Grid1.Node(i1, i2, i3, DV_ZERO_M1_ZERO) = Grid1.Node(i1, i2+1, i3, DV_ZERO_M1_ZERO);
        Grid1.Node(i1, i2, i3, DV_ZERO_ZERO_M1) = Grid1.Node(i1, i2, i3+1, DV_ZERO_ZERO_M1);


        Grid1.Node(i1, i2, i3, DV_M1_M1_ZERO) = Grid1.Node(i1+1, i2+1, i3, DV_M1_M1_ZERO);
        Grid1.Node(i1, i2, i3, DV_P1_M1_ZERO) = Grid1.Node(i1-1, i2+1, i3, DV_P1_M1_ZERO);
        Grid1.Node(i1, i2, i3, DV_ZERO_M1_M1) = Grid1.Node(i1, i2+1, i3+1, DV_ZERO_M1_M1);
        Grid1.Node(i1, i2, i3, DV_ZERO_P1_M1) = Grid1.Node(i1, i2-1, i3+1, DV_ZERO_P1_M1);
        Grid1.Node(i1, i2, i3, DV_M1_ZERO_M1) = Grid1.Node(i1+1, i2, i3+1, DV_M1_ZERO_M1);
        Grid1.Node(i1, i2, i3, DV_P1_ZERO_M1) = Grid1.Node(i1-1, i2, i3+1, DV_P1_ZERO_M1);

        Grid1.Node(i1, i2, i3, DV_M1_M1_M1) = Grid1.Node(i1+1, i2+1, i3+1, DV_M1_M1_M1);
        Grid1.Node(i1, i2, i3, DV_P1_M1_M1) = Grid1.Node(i1-1, i2+1, i3+1, DV_P1_M1_M1);
        Grid1.Node(i1, i2, i3, DV_P1_P1_M1) = Grid1.Node(i1-1, i2-1, i3+1, DV_P1_P1_M1);
        Grid1.Node(i1, i2, i3, DV_M1_P1_M1) = Grid1.Node(i1+1, i2-1, i3+1, DV_M1_P1_M1);
      }
}