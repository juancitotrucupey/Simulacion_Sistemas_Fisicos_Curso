#include<iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;


#define Lx 16
#define Nx 8 //# hilos por bloque

const int Mx=(Lx+Nx-1)/Nx; //# bloques por grilla

//__constant__ float d_w[5];

__global__ void SumaDosVectores(float *d_a, float *d_b,float *d_c){
int ix=blockIdx.x*blockDim.x+threadIdx.x
  d_c[ix]=d_a[ix]+d_b[ix];

}

int main(void){
  int ix;
  //Declarar las matrices
  //En el Host
  float h_a[Lx], h_b[Lx], h_c[Lx];
  //En el Device
  float *d_a; cudaMalloc((void**)&d_a,Lx*sizeof(float));
  float *d_b; cudaMalloc((void**)&d_b,Lx*sizeof(float));
  float *d_c; cudaMalloc((void**)&d_c,Lx*sizeof(float));
  /*
  //Enviarlas al device
  cudaMencpyToSymbol(d_w,h_w,5*sizeof(float),0,cudaMencpyHostToDevice);
  */
  //Inicializar los datos
  //cargar los datos en el host
  for(ix=0;ix<x;ix++)h_a[ix]=ix;
    for(ix=0;ix<x;ix++)h_b[ix]=2*ix;
    
  //Enviar al device
  cudaMencpy(d_a,h_a,Lx*sizeof(float),cudaMencpyHostToDevice);
  cudaMencpy(d_b,h_b,Lx*sizeof(float),cudaMencpyHostToDevice);


  //Procesar en la tarjeta grafica
   dim3 ThreadsPerBlock(Nx,1,1);  //x,y,z
  dim3 BlocksPerGrid(Mx,1,1);
  SumaDosVectores<<<BlocksPerGrid,ThreadsPerGrid>>>(d_a,d_b,d_c);

  //Devolver e imprimir
  //Devolverlos al Host
     cudaMencpy(d_c,h_c,Lx*sizeof(float),cudaMencpyDeviceToHost);
 
   //Imprimirlos
     for(ix=0;ix<Lx;ix++) cout<<h_c[ix]<<endl;


				       //Liberar la memoria				       
				       cudaFree(d_a); cudaFree(d_b); cudaFree(d_c);

				       return 0;

}
