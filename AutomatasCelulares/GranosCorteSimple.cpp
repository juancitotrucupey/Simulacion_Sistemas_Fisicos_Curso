// UnPlaneta por LeapFrog
#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;

const int L=10;
const double p=0.5;
class LatticeGas{
private:

int n[L][2], nnew[L][2];

public:
void Inicie(void);
void Muestre(void);
void Colision(void);
void 
};
void LatticeGas:: Inicie(void){
for(int ix=0;ix<L;ix++)
for(int i=0;i<L;i++)
n[ix][i]=0;
n[L/2][0]=1;
}
void LatticeGas::Muestre(void){
for(int i=0;i<2;i++){
for(int ix=0;ix<L;ix++)
cout<<n[ix][i]<<" ";
cout<<endl;
}
}

void LatticeGas::Colision(Crandom & ran){
for(int i=0;i<2;i++){
for(int ix=0;ix<L;ix++)

int main(void){
LatticeGas Difusion;
Crandom ran(12);
Difusion.Inicie();
Difusion.Muestre();
return 0;
}
