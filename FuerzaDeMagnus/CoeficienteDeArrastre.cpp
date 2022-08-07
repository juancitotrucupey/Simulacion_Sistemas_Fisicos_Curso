// Efecto Magnus
#include <iostream>
#include <cmath>
#include "Vector.h"
using namespace std;

int main(void){

double Re,Ca;

int i,j,N=4;

double XData [N];
double YData [N];
XData[0]=pow(10,5)*(pow(2,0.3)-pow(2,0.575)); YData[0]=0.22;
XData[1]=0; YData[1]=0;
XData[2]=pow(10,5)*(pow(2,0.85)-pow(2,0.575)); YData[2]=-0.22;
XData[3]=pow(10,5)*(pow(2,0.575)-pow(2,0.85)); YData[3]=0.22;
//Parametros del polinomio
double NRe;
double Pol=0;
double Mul=1;

//Parametros del ajuste exponencial

double A, Alpha,SRe;


for(Re=1000;Re<1000000;Re+=1){


if(Re < 100000*pow(2,0.3)){
Ca=0.5;
}

else if((100000*pow(2,0.3)<= Re)&(Re <=100000*pow(2,0.85))){
NRe=Re-pow(2,0.575)*pow(10,5);


for(i=0;i<4;i++){
Mul=Mul*YData[i];
for(j=0;j<4;j++){
if(i!=j){
Mul=Mul*(NRe-XData[j])/(XData[i]-XData[j]);
}
}
Pol+=Mul;
Mul=1;
}
Ca =Pol+0.28;
Pol=0;
}

else if(Re> 100000*pow(2,0.85)){

A=0.375-0.06;
Alpha=log(A/(0.375-0.372))/( pow(10,5)*(10-pow(2,0.85)));
SRe=Re-100000*pow(2,0.85);
Ca=0.375-A*exp((-1)*Alpha*SRe);

}

cout<<Re<<" "<<Ca<<endl;
}




  return 0;
}
