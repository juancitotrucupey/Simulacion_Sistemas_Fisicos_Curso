#include <iostream>
#include <cmath>

using namespace std;

double f (double x, double t)
{

  return x;
}
double Xexacta(double t)
{
  return exp(t);
}



void UnPasoRungeKutta(double & x, double & t, double dt)
{ double dx1, dx2, dx3, dx4;
  double dx;
  dx1=dt*f(x,t);
  dx2 = dt*f(x+dx1/2,t+dt/2);
  dx3 = dt*f(x+dx2/2,t+dt/2);
  dx4 = dt*f(x+dx3,t+dt);
  
  x+=(dx1 + 2*dx2 + 2*dx3 + dx4)/6;
  t+=dt;
}

int main(void){

  double t=0, x=1;
  double dt=0.01;
  
  for(t=0,x=1;t<5;)
    {
      cout<<t<<" "<<x<<" "<<Xexacta(t)<<endl;
      UnPasoRungeKutta(x,t,dt);
    }

  
  return 0;
}
