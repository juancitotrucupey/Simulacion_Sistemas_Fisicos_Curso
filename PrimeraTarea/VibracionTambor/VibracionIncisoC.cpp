#include <iostream>
#include <cmath>

using namespace std;



double f1 (double x1 , double x2, double r, double Lambda)
{
  return x2;
}
double f2(double x1, double x2,double r, double Lambda)
{
  return -(r*x2 + pow(Lambda,2)*pow(r,2)*x1)/pow(r,2);
}



void UnPasoRungeKutta4Acoplado(double & x1, double & x2, double & r, double dr, double Lambda)
{ double dx11, dx12, dx13, dx14;                      double dx21, dx22, dx23, dx24;
  
  dx11=dr*f1(x1,x2,r,Lambda);                                dx21=dr*f2(x1,x2,r,Lambda);
  dx12 = dr*f1(x1+dx11/2,x2+dx21/2,r+dr/2,Lambda);           dx22 = dr*f2(x1+dx11/2,x2+dx21/2,r+dr/2,Lambda);
  dx13 = dr*f1(x1+dx12/2,x2+dx22/2,r+dr/2,Lambda);           dx23 = dr*f2(x1+dx12/2,x2+dx22/2,r+dr/2,Lambda);
  dx14 = dr*f1(x1+dx13,x2 + dx23,r+dr,Lambda);               dx24 = dr*f2(x1+dx13,x2 + dx23,r+dr,Lambda);
  
  x1+=(dx11 + 2*dx12 + 2*dx13 + dx14)/6;               x2+=(dx21 + 2*dx22 + 2*dx23 + dx24)/6;
  r+=dr;
}

double FuncionRen1 (double Lambda)
{

  double r,x1,x2;
  double dr=0.001;
  
  for(r=0.01,x1=1,x2=0;r<=1.0;)
    {
      UnPasoRungeKutta4Acoplado(x1,x2,r,dr,Lambda);
    }
  

  return x1;
}

void CeroPorBiseccion(double & Lambda,double a, double b)
{
  double fa=FuncionRen1(a),fb=FuncionRen1(b);
  double m,fm;
  double Err=1e-7;
    
  while((b-a)>Err)
    {
      m=(a+b)/2;
      fm=FuncionRen1(m);
      if(fb*fm<=0)
	{a=m;
	  fa=fm;}
      else
	{b=m;
	  fb=fm;}
      
    }

  Lambda = (a+b)/2;

}

int main(void){
  double Lambda1,  Lambda2,  Lambda3,  Lambda4,  Lambda5;

  CeroPorBiseccion(Lambda1,2.4,2.45); CeroPorBiseccion(Lambda2,5,6); CeroPorBiseccion(Lambda3,8.4,8.8); CeroPorBiseccion(Lambda4,11.6,11.9); CeroPorBiseccion(Lambda5,14.9,14.97);
    cout<<Lambda1<<" "<<Lambda2<<" "<<Lambda3<<" "<<Lambda4<<" "<<Lambda5<<endl;

  double x11,x12,x21,x22,x31,x32,x41,x42,x51,x52;
  double r1,r2,r3,r4,r5;
  double dr=0.001;

 
  for(r1=0.0001,r2=0.0001,r3=0.0001,r4=0.0001,r5=0.0001,x11=1,x21=1,x31=1,x41=1,x51=1,x12=0,x22=0,x32=0,x42=0,x52=0;r1<=1.0;)
    {
      cout<<r1<<" "<<x11<<" "<<x21<<" "<<x31<<" "<<x41<<" "<<x51<<endl;
      
      UnPasoRungeKutta4Acoplado(x11,x12,r1,dr,Lambda1);
 UnPasoRungeKutta4Acoplado(x21,x22,r2,dr,Lambda2);
     UnPasoRungeKutta4Acoplado(x31,x32,r3,dr,Lambda3);
 UnPasoRungeKutta4Acoplado(x41,x42,r4,dr,Lambda4);
  UnPasoRungeKutta4Acoplado(x51,x52,r5,dr,Lambda5);
    }


 
  
  return 0;
}
