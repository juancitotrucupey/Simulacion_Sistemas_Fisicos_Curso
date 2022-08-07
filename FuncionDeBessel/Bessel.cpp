#include <iostream>
#include <cmath>

using namespace std;

double f(double t, double alpha, double x)
{ return cos(alpha*t - x*sin(t));
}

double IntegralPorSimpson(double alpha, double x, double a, double b, int N)
{
double h=(b-a)/N;
int i; double ti,suma,integral;

suma=0;
for(i=0;i<=N;i++)
    {
    ti=a+i*h;
    if((i==0)||(i==N))
      suma+=f(ti,alpha,x);
    else if(i%2==0)
      suma+=2*f(ti,alpha,x);
    else
      suma+=4*f(ti,alpha,x);
}
 integral=(h/3)*suma;
return integral;
}

double Bessel (double alpha, double x)
{
  return (1.0/M_PI)*IntegralPorSimpson(alpha,x,0,M_PI,100);
}
double CeroPorBiseccion(double a, double b,double alpha)
{
  double fa=Bessel(alpha,a),fb=Bessel(alpha,b);
  double m,fm;
  double Err=1e-7;
  
  while((b-a)>Err)
    {
      m=(a+b)/2;
      fm=Bessel(alpha,m);
      if(fb*fm<=0)
	{a=m;
	  fa=fm;}
      else
	{b=m;
	  fb=fm;}
      
    }

      return (a+b)/2;

}

int main(void)
{
  double alpha=4;
  double a=0;
  double b=5;

  cout<<CeroPorBiseccion(a,b,alpha)<<endl;
     
return 0;

}
