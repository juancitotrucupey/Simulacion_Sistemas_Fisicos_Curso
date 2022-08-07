// Colisionador de billar
#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;
//Constantes del integrador
const double CHI   =0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double XI    =-0.06626458266981849;
const double uno_m2LAMBDA_S2=(1-2*LAMBDA)/2;
const double uno_m2_XIplusCHI=(1-2*(XI+CHI));


//Constantes del programa
const double ERFF=1e-9;
const double Deltat=1e-4;
const double g=9.8;
const int N=5; //Cantidad de bolas interactuando
//_________________Posición de las bandas_______
//Norte
const double BN=0.50;
//Sur
const double BS=0;
// Este
  const double BE=0.50;
//Oeste
const double BW=0;


//Constantes de la simmulación de billar, tomado de http://billiards.colostate.edu/threads/physics.html 
const double ResBallTable = 0.5;
const double ResBallBall = 0.94;
const double BallRadious = 0.0571/2;
const double BallMass = 0.170097;
const double MuBallRoll = 0.01;
const double MuBallBall = 0.04;
const double MuBallTableSlid = 0.2;
const double UssualRadiousContact = 2e-3;
const double BallBallCollTime = 270e-6;
const double BBInterDist = 6e-4;
//Constantes físicas de la simulación
const double Inercia=2*BallMass*BallRadious*BallRadious/5;
//______Ball-Table Normal Interaction
const double BTInterDist = UssualRadiousContact*UssualRadiousContact/BallRadious;
const double BTElasticModuli = 3*BallMass*g/(4*UssualRadiousContact*BTInterDist);
const double BTNorElastick = 2*BTElasticModuli*UssualRadiousContact;
const double BTCollTime = M_PI*sqrt(BallMass/(2*BTNorElastick));
const double BTDamping = -2*log(ResBallTable)/BTCollTime;
//____________Ball-Ball Normal Interaction
const double BBNorElastick = M_PI*M_PI*BallMass/(2*BallBallCollTime*BallBallCollTime);
const double BBDamping = -2*log(ResBallBall)/BallBallCollTime;
const double BBElasticModuli = BBNorElastick/(2*sqrt(BallRadious*BBInterDist/2));
//____________Ball-Table Tangent Interaction
const double BTKundall = BTNorElastick;
const double BTTanDamping = BTDamping/2;
//___________Ball-Ball Tangent Interaction
const double BBKundall = BBNorElastick;
const double BBTanDamping = BBDamping/2;



class Cuerpo;
class Colisionador;


/*La informción de contacto está guardada en cada partícula, no pertenece al colisionador, y como el colisionador es una clase amiga puede acceder a est información y modificarla, lo cual es lo que se hacia con la disposición inicial del programa*/
//------------------ Clase Cuerpo -----------------------
class Cuerpo{
private:
  vector3D r,V,Ftot,Fnor,Ftan,Ffloor,Omega,Tau,dcontacto[N+5];
  double m,R,I,hold[N+5];
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0
	      ,double omegax0,double omegay0,double omegaz0,double m0,double R0);
  void BorreFuerzaYTorque(void);
  void AgregueFuerzaGravedad(void);
  void AdicioneFuerza(vector3D F0);
  void AdicioneTorque(vector3D Tau0);
  void Mueva_r(double dt,double CTE);
  void Mueva_V(double dt,double CTE);
  void Dibujese(void);
  void RenueveContacto(vector3D dcontacto1,double hold1, int k);
  double Getx(void){return r.x();}; //Funcion inline (macro)
  double Gety(void){return r.y();}; //Funcion inline (macro)
  double Getz(void){return r.z();}; //Funcion inline (macro)
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0
		    ,double omegax0,double omegay0,double omegaz0, double R0, double m0){
  r.cargue(x0,y0,z0);   V.cargue(Vx0,Vy0,Vz0);    Omega.cargue(omegax0,omegay0,omegaz0);
  m=m0;  R=R0;   I=2*m0*R0*R0/5;
  for(int i=0;i<N+5;i++){
    dcontacto[i].cargue(0,0,0);
    hold[i]=0;
  }
}
void Cuerpo::BorreFuerzaYTorque(void){
  Ftot.cargue(0,0,0);   Tau.cargue(0,0,0); 
}
void Cuerpo::AgregueFuerzaGravedad(void){
  vector3D Fg; Fg.cargue(0,0,-m*g); Ftot+=Fg;
}
void Cuerpo::AdicioneFuerza(vector3D F0){
  Ftot+=F0;
}
void Cuerpo::AdicioneTorque(vector3D Tau0){
  Tau+=Tau0;
}
void Cuerpo::Mueva_r(double dt,double CTE){
  r+=V*(CTE*dt); 
}
void Cuerpo::Mueva_V(double dt,double CTE){
  V+=Ftot*(CTE*dt/m);  Omega+=Tau*(CTE*dt/I);
}
void Cuerpo::Dibujese(void){
 cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
void Cuerpo::RenueveContacto(vector3D dcontacto1,double hold1,int k){
  dcontacto[k]=dcontacto1;
  hold[k]=hold1;
}

//------------------ Clase Colisionador -----------------------
class Colisionador{
private:
 
 public:
  void CalculeTodasLasFuerzas(Cuerpo * Grano,double dt,double t);
  void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2,int j,double dt);
  void CalculeFuerzaFrontera(Cuerpo & Grano,int i,double dt);
 };

void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Grano,double dt, double t){
  int i,j,k;
  vector3D dcontactoBall;
  //Borrar las fuerzas de todos los planetas los Granos
  for(i=0;i<N;i++)    
    Grano[i].BorreFuerzaYTorque();
  //Agregar Fuerza de Gravedad
  for(i=0;i<N;i++)    
    Grano[i].AgregueFuerzaGravedad();
 
  //Hacer todas las bolas y calcular la fuerza entre estas
 
  for(i=0;i<N;i++){
       for(j=i+1;j<N;j++){
	 dcontactoBall=Grano[i].dcontacto[j];
	 cout<<"TIEMPO "<<t<<endl;
	 cout<<"Distancia de contacto entre bolas "<<i<<", "<<j<<" antes de calcular fuerza"<<endl;
	 cout<<dcontactoBall.x()<<" "<<dcontactoBall.y()<<" "<<dcontactoBall.z()<<endl;
	 CalculeFuerzaEntre(Grano[i],Grano[j],j,dt);
	  dcontactoBall=Grano[i].dcontacto[j];
	  cout<<"Distancia de contacto entre bolas "<<i<<", "<<j<<" despues de calcular fuerza"<<endl;
	 cout<<dcontactoBall.x()<<" "<<dcontactoBall.y()<<" "<<dcontactoBall.z()<<endl;
       }
  }
   if(t==0) cout<<"orden de las fronteras, z, BE,BW,BN,BS"<<endl;
  //Calcular la fuerza entre las bolas y las fronteras
  for(i=0;i<N;i++){
       for(k=0;k<5;k++){
	  dcontactoBall=Grano[i].dcontacto[N+k];
	 cout<<"Distancia de contacto entre bola "<<i<<" y las fronteras antes de calcular fuerza"<<endl;
		 cout<<dcontactoBall.x()<<" "<<dcontactoBall.y()<<" "<<dcontactoBall.z()<<endl; 
CalculeFuerzaFrontera(Grano[i],k,dt);
  dcontactoBall=Grano[i].dcontacto[N+k];
   cout<<"Distancia de contacto entre bola "<<i<<" y las fronteras despues de calcular fuerza"<<endl;
		 cout<<dcontactoBall.x()<<" "<<dcontactoBall.y()<<" "<<dcontactoBall.z()<<endl; 
	     }
    
  }
  }

  void Colisionador::CalculeFuerzaFrontera(Cuerpo & Grano,int i,double dt){
    
    vector3D r,F2,Fn,Ft,Vc,Vcn,Vct,n,t,t1,Vel,Omega;
  vector3D Trolling_resistance,Frolling_resistance,dcontactoFrontera;
  double d,h,normaFn,Ftmax,m1,R1,normaVcn,normaVct,normaFt,normaV,normaOmega,holdFrontera;
  dcontactoFrontera=Grano.dcontacto[N+i];
  holdFrontera=Grano.hold[N+i];

  if(i==0){
    // cout<<"z"<<endl;
     r.cargue(0,0,Grano.Getz());
     h=Grano.R-Grano.Getz();
  }
  else if(i==1){
    // cout<<"Be"<<endl;
     r.cargue(Grano.Getx()-BE,0,0);
      h=Grano.R+Grano.Getx()-BE;
  }
  else if(i==2){
    // cout<<"BW"<<endl;
    r.cargue(Grano.Getx()-BW,0,0);
     h=BW-(Grano.Getx()-Grano.R);
  }
  else if(i==3){
    //cout<<"BN"<<endl;
    r.cargue(0,Grano.Gety()-BN,0);
     h=Grano.R+Grano.Gety()-BN;
  }
  else if(i==4){
    //cout<<"BS"<<endl;
    r.cargue(0,Grano.Gety()-BS,0);
     h=BS-(Grano.Gety()-Grano.R);
  }
 
      d=norma(r);
      
      if(h>0){
      
    //Geometría y dinámica del contacto
    m1=Grano.m;  
    R1=Grano.R;
    n=r/d;
    //Calcular velocidad de contacto y el vector tangente
    Vc=Grano.V+(Grano.Omega^n)*(-R1);
    normaVcn=Vc*n; Vcn=n*normaVcn; Vct=Vc-Vcn;  normaVct=norma(Vct);

    //FUERZAS NORMALES
    //Fuerza de Hertz
   
     normaFn=4*BTElasticModuli*sqrt(R1)*pow(h,1.5)/3; 
    //Disipacion
    //__________________________________________________________________
    //Damping independent of interpenetration
    normaFn-= m1*BTDamping*normaVcn; if(normaFn<0) normaFn=0;
    //__________________________________________________________________
    //Damping dependent of interpenetration
    // normaFn-= m1*sqrt(h)*Gamma*normaVcn; if(normaFn<0) normaFn=0;
    //__________________________________________________________________   

    Fn=n*normaFn;

    //Rolling resistance
     if(normaVct<ERFF){
      t.cargue(0,0,0);
      normaV=norma(Grano.V);
      if(normaV<ERFF) Vel.cargue(0,0,0); else Vel=Grano.V/normaV;
      Frolling_resistance = (-normaFn)*MuBallRoll*Vel;	
      Grano.AdicioneFuerza(Frolling_resistance);
      normaOmega=norma(Grano.Omega);
      if(normaOmega<ERFF){
	Trolling_resistance.cargue(0,0,0); 
	Grano.AdicioneTorque(Trolling_resistance);}
      else{ Omega=Grano.Omega/normaOmega;
      
      vector3D Uno = (Grano.Omega^n)*(-R1);
      double Uno_1 = norma(Uno);
       double Refectivo = Uno_1/normaOmega;
      
      Trolling_resistance = (-normaFn)*MuBallRoll*Inercia/(m1*Refectivo)*Omega; 
	Grano.AdicioneTorque(Trolling_resistance);
	 dcontactoFrontera.cargue(0,0,0);
      }
    }
    else{ t=Vct/normaVct;
      
    //FUERZAS TANGENCIALES
    // Resorte de kundall
      /*
dcontactoFrontera+=(Vct*dt);
    Ft=dcontactoFrontera*(-BTKundall);
    //fuerza cinética
    Ftmax=MuBallTableSlid*normaFn; normaFt=norma(Ft);
    if(normaFt>Ftmax) Ft=t*(-Ftmax);
      */       
     Ft=dcontactoFrontera*(-BBKundall);
    
    Ftmax=MuBallTableSlid*normaFn; normaFt=norma(Ft);

    if(normaFt<Ftmax){
      dcontactoFrontera+=(Vct*dt);
    }

    Ft=dcontactoFrontera*(-BBKundall);

    normaFt=norma(Ft);

    if(normaFt<ERFF) t1.cargue(0,0,0); else t1=Ft/normaFt;

    if(normaFt>=Ftmax){
      if(normaVct<ERFF) Ft=t1*(Ftmax); else Ft=(-Ftmax)*t;
      } 
      }
     //Torque sobre el eje de contacto
     double torque,Alpha;
      torque=MuBallTableSlid*normaFn*2*sqrt(BallRadious*h)/3;
      normaOmega=norma(Grano.Omega);
      if(normaOmega<ERFF){
	Omega.cargue(0,0,0);
	Alpha=0;} 
      else{ Omega=Grano.Omega/normaOmega;
	if(h<ERFF)Alpha=0;
	else Alpha=log(2)/(2*BallRadious*h*normaOmega);}
      double coseno; coseno=n*Omega;
      double FuncionVel; FuncionVel=exp((-1)*Alpha*pow(normaVct,2));
      vector3D Torque; Torque=(-1)*coseno*torque*FuncionVel*n;
       Grano.AdicioneTorque(Torque);
    
  //Construir la fuerza total
      
F2=Fn+Ft;
Grano.AdicioneFuerza(F2);
 Grano.AdicioneTorque((n*(-R1))^Ft);      
   
  }
  else if(holdFrontera>=0){
    //Reiniciar dcontacto en ceros
    dcontactoFrontera.cargue(0,0,0);
  }
      
      Grano.RenueveContacto(dcontactoFrontera,holdFrontera,N+i);
    }

void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2,int j,double dt){
  vector3D r21,F2,Fn,Ft,Vc,Vcn,Vct,n,t,t1,dcontactoBall,Omega; 
  double d21,h,normaFn,Ftmax,m1,m2,m12,R1,R2,normaVcn,normaVct,normaFt,hold,normaOmega;
  dcontactoBall=Grano1.dcontacto[j];
  hold=Grano1.hold[j];
  r21=Grano2.r-Grano1.r; d21=norma(r21); h=(Grano2.R+Grano1.R)-d21;

  if(h>0){
    //Geometría y dinámica del contacto
    m1=Grano1.m;   m2=Grano2.m;   m12=(m1*m2)/(m1+m2);
    R1=Grano1.R;   R2=Grano2.R;
   n=r21/d21;
    //Calcular velocidad de contacto y el vector tangente
    Vc=(Grano2.V-Grano1.V)+(((Grano2.Omega*(-R2))-(Grano1.Omega*R1))^n);
    normaVcn=Vc*n; Vcn=n*normaVcn; Vct=Vc-Vcn;  normaVct=norma(Vct);
    if(normaVct<ERFF) t.cargue(0,0,0); else t=Vct/normaVct;

    //FUERZAS NORMALES
    //Fuerza de Hertz
    double RadioEfectivo=Grano2.R*Grano1.R/(Grano2.R+Grano1.R);
     normaFn=4*BBElasticModuli*sqrt(RadioEfectivo/2)*pow(h,1.5)/3; 
    
    //Disipacion
    //__________________________________________________________________
    //Damping independent of interpenetration
    normaFn-= m12*BBDamping*normaVcn; if(normaFn<0) normaFn=0;
    //__________________________________________________________________
    //Damping dependent of interpenetration
    // normaFn-= m12*sqrt(h)*Gamma*normaVcn; if(normaFn<0) normaFn=0;
    //__________________________________________________________________   
    
    Fn=n*normaFn;
    
    //FUERZAS TANGENCIALES
    
    // Resorte de kundall
 //fuerza estática
    /* 
    dcontactoBall+=(Vct*dt);
    Ft=dcontactoBall*(-BBKundall);
    //fuerza cinética
    Ftmax=MuBallBall*normaFn; normaFt=norma(Ft);
    if(normaFt>Ftmax) Ft=t*(-Ftmax);
    */
    //Si se cambia el  

    Ft=dcontactoBall*(-BBKundall);
    Ftmax=MuBallTableSlid*normaFn; normaFt=norma(Ft);
   
    if(normaFt<Ftmax){
      dcontactoBall+=(Vct*dt);
    }

    Ft=dcontactoBall*(-BBKundall);

    normaFt=norma(Ft);
    
    if(normaFt<ERFF) t1.cargue(0,0,0); else t1=Ft/normaFt;

    if(normaFt>=Ftmax){
      if(normaVct<ERFF) Ft=t1*(Ftmax); else Ft=t*(-Ftmax);
      }
     //Torque sobre el eje de contacto
    double torque,Alpha;
    vector3D Omega1,Omega2;
    Omega1=(n*Grano1.Omega)*n; Omega2=(n*Grano2.Omega)*n;
      Omega=Omega2-Omega1;
      normaOmega=norma(Omega);
      if(normaOmega<ERFF){
	Omega.cargue(0,0,0);
	Alpha=0;} 
      else{ Omega=Omega/normaOmega;
	if(h<ERFF)Alpha=0;
	else Alpha=log(2)/(2*RadioEfectivo*h*normaOmega);}
      double coseno; coseno=n*Omega;
      double FuncionVel; FuncionVel=exp((-1)*Alpha*pow(normaVct,2));
      vector3D Torque; Torque=(-1)*coseno*torque*FuncionVel*n;
       Grano2.AdicioneTorque(Torque);
       Grano1.AdicioneTorque((-1)*Torque);
//Construir la fuerza total
        
F2=Fn+Ft;
    Grano2.AdicioneFuerza(F2);
    Grano2.AdicioneTorque((n*(-R2))^Ft);      
    Grano1.AdicioneFuerza(F2*(-1));
    Grano1.AdicioneTorque((n*R1)^(Ft*(-1)));      


  }
  else if(hold>=0){
    //Reiniciar dcontacto en ceros
    dcontactoBall.cargue(0,0,0);
  }
  hold=h;

  Grano1.RenueveContacto(dcontactoBall,hold,j);
}


//-------------------------------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'Billar2D.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-0.1:0.6]"<<endl;
  cout<<"set yrange [-0.1:0.6]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<0.5/7<<"*t,0";
    cout<<" , "<<0.5/7<<"*t,0.5";
    cout<<" , 0,"<<0.5/7<<"*t";
    cout<<" , 0.5,"<<0.5/7<<"*t";
}
void TermineCuadro(void){
    cout<<endl;
}



int main(void){
  Cuerpo Grano[N]; 
  Colisionador Hertz;
     double t;
  double tmax=1, tdibujo; int ndibujos=10000;
  double Deltat = 0.0001;
  int i;
 
  //Forma en que se ponen las condiciones iniciales
//Grano[0].Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double omegax0,double omegay0,double omegaz0, double R0, double m0);  

  Grano[0].Inicie(0.1,0.2,BallRadious,0,0,-2,100,0,0,BallRadious,BallMass);
  Grano[1].Inicie(0.4,0.2,BallRadious,10,0,0,0,0,0,BallRadious,BallMass);
   Grano[2].Inicie(0.25,0.2,BallRadious,-3,-2,0,1,100,1,BallRadious,BallMass);
  Grano[3].Inicie(0.4,0.1,BallRadious,-3,-2,0,100,1,1,BallRadious,BallMass);
  Grano[4].Inicie(0.25,0.4,BallRadious,-3,-2,0,10,10,10,BallRadious,BallMass);
  
   InicieAnimacion();
    
  for(t=0,tdibujo=0;t<tmax;t+=Deltat,tdibujo+=Deltat){
    //Animacion
    /*
    if(tdibujo>tmax/ndibujos){
      InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
      }*/
    /*
    InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
    */
    for(i=0;i<=N;i++) Grano[i].Mueva_r(Deltat,CHI);
    Hertz.CalculeTodasLasFuerzas(Grano,Deltat,t); for(i=0;i<=N;i++) Grano[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
    for(i=0;i<=N;i++) Grano[i].Mueva_r(Deltat,XI);
    Hertz.CalculeTodasLasFuerzas(Grano,Deltat,t); for(i=0;i<=N;i++) Grano[i].Mueva_V(Deltat,LAMBDA); 
    for(i=0;i<=N;i++) Grano[i].Mueva_r(Deltat,uno_m2_XIplusCHI);
    Hertz.CalculeTodasLasFuerzas(Grano,Deltat,t); for(i=0;i<=N;i++) Grano[i].Mueva_V(Deltat,LAMBDA); 
    for(i=0;i<=N;i++) Grano[i].Mueva_r(Deltat,XI);
    Hertz.CalculeTodasLasFuerzas(Grano,Deltat,t); for(i=0;i<=N;i++) Grano[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
    for(i=0;i<=N;i++) Grano[i].Mueva_r(Deltat,CHI);

    
  }

  return 0;
}
