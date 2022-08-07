#include <iostream>
#include <cmath>

using namespace std;

double f(double a)
{ return sin(a);
}

double IntegralPorSimpson( double a, double b, int N)
{
double h=(b-a)/N;
int i; double xi,suma,integral;

suma=0;
for(i=0;i<=N;i++)
    {
    xi=a+i*h;
    if((i==0)||(i==N))
    suma+=f(xi);
    else if(i%2==0)
    suma+=2*f(xi);
    else
suma+=4*f(xi);
}
integral=h/3*suma;
return integral;
}


int main(void)
{
double a=0,b=M_PI; int N=20;
double h=(b-a)/N;

cout<<h<<" "<<IntegralPorSimpson(a,b,N)<<endl;

return 0;

}
