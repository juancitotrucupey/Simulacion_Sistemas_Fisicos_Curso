#include<iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

#define Lx 16
#define Nx 8 //# de hilos por bloque
const int Mx=(Lx+Nx-1)/Nx; //# de bloques por grilla

//__constant__ float d_w[5];
__global__ void SumaDosVectores(float *d_a,float *d_b,float *d_c){
  int ix=blockIdx.x*blockDim.x+threadIdx.x;
  d_c[ix]=d_a[ix]+d_b[ix];
}

int main(void){
  int ix;
  //DECLARAR LAS MATRICES
  //En el Host
  float h_a[Lx],h_b[Lx],h_c[Lx];
  //En el Device
  float *d_a;  cudaMalloc((void**) &d_a,Lx*sizeof(float));
  float *d_b;  cudaMalloc((void**) &d_b,Lx*sizeof(float));
  float *d_c;  cudaMalloc((void**) &d_c,Lx*sizeof(float));

  /*
  //Cargar las constantes en el Host
  float h_w[5];
  h_w[0]=1.0/3;    h_w[1]=h_w[2]=h_w[3]=h_w[4]=1.0/6;
  //Enviarlas al Device
  cudaMemcpyToSymbol(d_w,h_w,5*sizeof(float),0,cudaMemcpyHostToDevice);
  */
  
  //INICIALIZAR LOS DATOS
  //Cargarlos en el Host
  for(ix=0;ix<Lx;ix++) h_a[ix]=ix;
  for(ix=0;ix<Lx;ix++) h_b[ix]=2*ix;
  //Enviarlos al Device
  cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_b,h_b,Lx*sizeof(float),cudaMemcpyHostToDevice);

  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,1,1); //x,y,z
  dim3 BlocksPerGrid(Mx,1,1);
  SumaDosVectores<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,d_b,d_c);

  //IMPRIMIR LOS DATOS
  //Devolverlos al Host
  cudaMemcpy(h_c,d_c,Lx*sizeof(float),cudaMemcpyDeviceToHost);
  //Imprimirlos
  for(ix=0;ix<Lx;ix++) cout<<h_c[ix]<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_a);  cudaFree(d_b);  cudaFree(d_c);

  return 0;
}
