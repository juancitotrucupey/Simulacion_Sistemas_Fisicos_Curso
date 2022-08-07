// UnPlaneta por LeapFrog
#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;

int main(void){
Crandom ran1(0), ran2(10);
int i;
double tau=3.0;
for(i=0;i<10000;i++){
cout<<sqrt(-2*log(ran1.r()))*cos(2*M_PI*ran2.r())<<endl;
}  
  return 0;
}
