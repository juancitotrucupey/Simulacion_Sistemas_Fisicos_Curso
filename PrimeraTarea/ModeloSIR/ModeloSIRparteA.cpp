#include <iostream>
#include <cmath>

using namespace std;

  
const double Omega0=1.0;

double f1 (double x1, double x2, double t,double b,double g)
{
  return -b*x1*x2;
}
double f2 (double x1, double x2,double t, double b, double g)
{
  return b*x1*x2 - g*x2;
}



void UnPasoRungeKutta4Acoplado(double & x1, double & x2, double & t, double dt, double b, double g)
{ double dx11, dx12, dx13, dx14;                      double dx21, dx22, dx23, dx24;
  
  dx11=dt*f1(x1,x2,t,b,g);                                dx21=dt*f2(x1,x2,t,b,g);
  dx12 = dt*f1(x1+dx11/2,x2+dx21/2,t+dt/2,b,g);           dx22 = dt*f2(x1+dx11/2,x2+dx21/2,t+dt/2,b,g);
  dx13 = dt*f1(x1+dx12/2,x2+dx22/2,t+dt/2,b,g);           dx23 = dt*f2(x1+dx12/2,x2+dx22/2,t+dt/2,b,g);
  dx14 = dt*f1(x1+dx13,x2 + dx23,t+dt,b,g);               dx24 = dt*f2(x1+dx13,x2 + dx23,t+dt,b,g);
  
  x1+=(dx11 + 2*dx12 + 2*dx13 + dx14)/6;               x2+=(dx21 + 2*dx22 + 2*dx23 + dx24)/6;
  t+=dt;
}

int main(void){

  double t,x1,x2,x3,b,g;
  double dt=0.1;
  b=0.24;
  g=0.08;
  
  for(t=0,x1=0.999,x2=0.001;f1(x1,x2,t,b,g)<-0.00001 || x1>0.8;)
    {
      x3=1-x2-x1;
      cout<<t<<" "<<x1<<" "<<x2<<" "<<x3<<" "<<f1(x1,x2,t,b,g)<<endl;
      UnPasoRungeKutta4Acoplado(x1,x2,t,dt,b,g);  
    }

  
  return 0;
}
