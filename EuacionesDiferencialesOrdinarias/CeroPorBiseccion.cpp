#include <iostream>
#include <cmath>

using namespace std;

const double D=10;
const double d=1;
const double pi=3.14159265359;


double dispersion(double t){
  double a,b,c;
  a=1/sqrt(2*D*d*t*2*pi);
  b=(-0.5)*pow(60/sqrt(2*D*d*t),2);
  c=a*exp(b)-pow(10,-3);
  return c;
}

double CeroPorBiseccion(double a, double b)
{
double fa=dispersion(a),fb=dispersion(b);
  double m,fm;
  double Err=1e-7;
  
  while((b-a)>Err)
    {
      m=(a+b)/2;
      fm=dispersion(m);
      if(fb*fm<=0)
	{a=m;
	  fa=fm;}
      else
	{b=m;
	  fb=fm;}
      
    }

      return (a+b)/2;

}

int main(void){

double a,b,c;
/*
cin>>a;
cin>>b;
cout<<dispersion(a)<<endl;
cout<<dispersion(b)<<endl;
if(dispersion(a)*dispersion(b)>=0){cin>>a; cin>>b;}

c=CeroPorBiseccion(a,b);
cout<<c<<endl;
*/
for(a=0.0000001;a<10000;a+=0.001){
cout<<a<<" "<<dispersion(a)<<endl;
}
return 0;
}
