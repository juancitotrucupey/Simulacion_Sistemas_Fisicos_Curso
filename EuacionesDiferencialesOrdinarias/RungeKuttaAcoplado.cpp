#include <iostream>
#include <cmath>

using namespace std;

const double Omega0=1.0;

double f1 (double x1, double x2, double t)
{
  return x2;
}
double f2(double x1, double x2,double t)
{
  return -Omega0*Omega0*x1;
}



void UnPasoRungeKutta4Acoplado(double & x1, double & x2, double & t, double dt)
{ double dx11, dx12, dx13, dx14;                      double dx21, dx22, dx23, dx24;
  
  dx11=dt*f1(x1,x2,t);                                dx21=dt*f2(x1,x2,t);
  dx12 = dt*f1(x1+dx11/2,x2+dx21/2,t+dt/2);           dx22 = dt*f2(x1+dx11/2,x2+dx21/2,t+dt/2);
  dx13 = dt*f1(x1+dx12/2,x2+dx22/2,t+dt/2);           dx23 = dt*f2(x1+dx12/2,x2+dx22/2,t+dt/2);
  dx14 = dt*f1(x1+dx13,x2 + dx23,t+dt);               dx24 = dt*f2(x1+dx13,x2 + dx23,t+dt);
  
  x1+=(dx11 + 2*dx12 + 2*dx13 + dx14)/6;               x2+=(dx21 + 2*dx22 + 2*dx23 + dx24)/6;
  t+=dt;
}

int main(void){

  double t,x1,x2;
  double dt=0.1;
  
  for(t=0,x1=0,x2=1;t<10;)
    {
      cout<<t<<" "<<x1<<" "<<endl;
      UnPasoRungeKutta4Acoplado(x1,x2,t,dt);
    }

  
  return 0;
}
