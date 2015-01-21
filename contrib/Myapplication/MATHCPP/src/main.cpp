
#include "MATH.hxx"
#include <stdio.h>
#include <cmath>
int main(){
  
  
  MATH a;
  double x=1.;
  double y=2.;
  double t=a.plus(x,y);
  double s=a.plus(t,sin(x));
  printf(" %f", s);
  
  
  
  
  
  
  
 return 0; 
}
