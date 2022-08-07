#include<iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

#define Lx 16
#define Ly 8
#define Nx 8
#define Ny 8
const int Mx=(Lx+Nx-1)/Nx;
const int My=(Ly+Ny-1)/Ny;


__device__ float SumeleUno(float x){
  return x+1;
}
__global__ void SumarConstante(float * d_a,size_t pitcha){
  int ix,iy; float *aux;
  ix=blockIdx.x*blockDim.x+threadIdx.x;  iy=blockIdx.y*blockDim.y+threadIdx.y;
  
  aux=d_a+(ix*pitcha)/sizeof(float)+iy; // aux es &(d_a[ix][iy])
  
  (*aux)=SumeleUno((*aux)); //  (*aux) es d_a[ix][iy]
}

//---------------- Programa Principal y Funciones Globales -------
int main(void){
  int ix,iy;



  //DECLARAR LAS MATRICES
  //En el Host
  float h_a[Lx][Ly];
  //En el Device
  float*d_a; size_t pitcha; cudaMallocPitch((void**) &d_a,&pitcha,Ly*sizeof(float),Lx);

  //INICIALIZAR LOS DATOS
  //Cargar los datos en el Host
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      h_a[ix][iy]=Ly*ix+iy;
  //Mostrar
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++)
      cout<<h_a[ix][iy]<<" ";
    cout<<endl;
  }
  cout<<endl;
  //Enviarlos al Device
  cudaMemcpy2D(d_a,pitcha,h_a,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);

  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,Ny,1); //x,y,z
  dim3 BlocksPerGrid(Mx,My,1);
  SumarConstante<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,pitcha);

  //IMPRIMIR LOS DATOS
  //Devolverlos al Host
  cudaMemcpy2D(h_a,Ly*sizeof(float),d_a,pitcha,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  //Imprimirlos
  //Mostrar
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++)
      cout<<h_a[ix][iy]<<" ";
    cout<<endl;
  }
  cout<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_a);

  return 0;
}
