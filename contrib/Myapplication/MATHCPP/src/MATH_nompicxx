#include "MATH.hxx"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
const double MATH::pi=3.1415926535897931;
const double MATH::epsilon=1.0E-12;



// ----------------------------------------------
MATH::MATH(){
//     int argc = 0;
//     char ** argv = NULL;
//     MPI_Init(&argc, &argv);
  return;
}
// ----------------------------------------------
MATH::~MATH() {
// MPI_Finalize();
return;
}



double MATH::sum(const std::vector<double>& _tab){
    double S=0.0;
    for (int i=0; i!=_tab.size(); ++i)	S+=_tab[i];
    return S;
}

double MATH::squareroot(double _x){
    return sqrt(_x);
}
double MATH::plus(double _x, double _y){
    return _x+_y;
}

double MATH::minus(double _x, double _y){
    return _x-_y;
}

double MATH::times(double _x, double _y){
    return _x*_y;
}
double MATH::sinx(double _x){
        return sin(_x);
}
