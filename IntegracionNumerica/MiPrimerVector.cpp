// Metodo de la lanzadera
#include <iostream>
#include <cmath>
#include "Vector.h"
using namespace std;

int main()
{
  vector3D a,b;
  a.cargue(1,2,3);
  a.show();

  b.cargue(5,6,7);
  cout<<a.operator*(b)<<endl;

  return 0;
  

}

// g++ UnPlaneta.cpp
// ./a.out | gnuplot
