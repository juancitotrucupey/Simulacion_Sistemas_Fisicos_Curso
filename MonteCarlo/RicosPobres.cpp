// UnPlaneta por LeapFrog
#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;

const int N=2000;


class Agente;
class Comercio;

//------------------ Clase Agente -----------------------
class Agente{
private:
  
  double capital;
public:
  void Inicie(double capital0);
  void AgregueGanancia(double ganancia);
  double GetCapital(void){return capital;}; //Funcion inline (macro)
  double GetApuesta(void); 
  friend class Comercio;
};
void Agente::Inicie(double capital0){
  capital = capital0;
}
void Agente::AgregueGanancia(double ganancia){
  capital+=ganancia;
}
double Agente::GetApuesta(void){
return 1.0;
}

//------------------ Clase Comercio -----------------------
class Comercio{
private:
  
public:
  void HacerUnaTransaccion(Agente * Empresa,Crandom & ran64);
  
};
 
void Comercio::HacerUnaTransaccion(Agente * Empresa, Crandom & ran64){
int i,j;
double bet_i,bet_j,toma_i,toma_j,ganancia_i,ganancia_j,nuevocapital_i,nuevocapital_j;

i=(int)N*ran64.r(); j=(int)(N-1)*ran64.r(); if(j>=i) j++; //Escoger dos agentes al azar;
bet_i=Empresa[i].GetApuesta(); bet_j=Empresa[j].GetApuesta(); //RecogerApuestas;
//CalcularCuantoGanaCadaUno
toma_i=(bet_i+bet_j)*ran64.r(); toma_j=(bet_i+bet_j)-toma_i;
ganancia_i=toma_i-bet_i; ganancia_j=toma_j-bet_j;
//RedistribuirElDinero
nuevocapital_i=Empresa[i].GetCapital()+ganancia_i;
nuevocapital_j= Empresa[j].GetCapital()+ganancia_j;
if(nuevocapital_i>0 & nuevocapital_j>0){
Empresa[i].AgregueGanancia(ganancia_i);
Empresa[j].AgregueGanancia(ganancia_j);
}

}




int main(void){
  Agente Empresa[N];
  Comercio BolsaDeValores;
  Crandom ran(0);
  int i,t,mcs,tmax=1000;
  double capital0=10;

//INICIE
for(i=0;i<N;i++){
Empresa[i].Inicie(10);}
//Hacer las transacciones
for(t=0;t<tmax;t++){
for(mcs=0;mcs<N;mcs++){
BolsaDeValores.HacerUnaTransaccion(Empresa,ran);
}}
//MUESTRE
for(i=0;i<N;i++) cout<<Empresa[i].GetCapital()<<endl;
  return 0;
}
