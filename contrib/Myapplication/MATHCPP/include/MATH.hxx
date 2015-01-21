#ifndef _MATH_HXX_
#define _MATH_HXX_

#include <vector>

class MATH
{
// Simple class
public:
   MATH();
  ~MATH();
    double sum(const std::vector<double>& tab);
    double squareroot(double x);
    double plus(double x, double y);
    double minus(double x, double y);
    double times(double x, double y);
    double sinx(double x);
private:
    static const double epsilon;
    static const double pi;
};

#endif
