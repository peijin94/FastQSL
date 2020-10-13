#include <math.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265
#endif

__device__ float get_Idx3d(float *Arr,long *AShapeN,long xIdx,long yIdx,long zIdx){
    return Arr[xIdx* AShapeN[1]*AShapeN[2]  +  yIdx* AShapeN[2]  +  zIdx];
}

__global__ void test_Idx2(float *Arr,long *AShapeN,long *getIdx,float *res){
    res[0] = get_Idx3d(Arr,AShapeN,getIdx[0],getIdx[1],getIdx[2]);
}



__global__ void Interp3d(float *Arr,long *AShapeN, float *inPoint,float *aget){

    //algorithm [https://core.ac.uk/download/pdf/44386053.pdf]

    float rx,ry,rz; // ratio of the point
    float Arr000,Arr001,Arr010,Arr011,Arr100,Arr101,Arr110,Arr111;
    float Aget;
    int x_Idx,y_Idx,z_Idx;

    rx = inPoint[0]-floorf(inPoint[0]);
    ry = inPoint[1]-floorf(inPoint[1]);
    rz = inPoint[2]-floorf(inPoint[2]);

    x_Idx = __float2int_rd(inPoint[0]);
    y_Idx = __float2int_rd(inPoint[1]);
    z_Idx = __float2int_rd(inPoint[2]);


    printf("%f %f %f\n", rx, ry,rz);
    printf("%f %f %f\n", inPoint[0], inPoint[1],inPoint[2]);
    printf("%d %d %d\n", x_Idx, y_Idx, z_Idx);


    Arr000 = get_Idx3d(Arr,AShapeN, x_Idx  ,y_Idx  ,z_Idx  );
    Arr001 = get_Idx3d(Arr,AShapeN, x_Idx  ,y_Idx  ,z_Idx+1);
    Arr010 = get_Idx3d(Arr,AShapeN, x_Idx  ,y_Idx+1,z_Idx  );
    Arr011 = get_Idx3d(Arr,AShapeN, x_Idx  ,y_Idx+1,z_Idx+1);
    Arr100 = get_Idx3d(Arr,AShapeN, x_Idx+1,y_Idx  ,z_Idx  );
    Arr101 = get_Idx3d(Arr,AShapeN, x_Idx+1,y_Idx  ,z_Idx+1);
    Arr110 = get_Idx3d(Arr,AShapeN, x_Idx+1,y_Idx+1,z_Idx  );
    Arr111 = get_Idx3d(Arr,AShapeN, x_Idx+1,y_Idx+1,z_Idx+1);
    //Arr000 = [x_Idx];

    Aget =      Arr000* (1.-rx)*(1.-ry)*(1.-rz)+\
                Arr001* (1.-rx)*(1.-ry)*(  rz)+\
                Arr010* (1.-rx)*(  ry)*(1.-rz)+\
                Arr011* (1.-rx)*(  ry)*(  rz)+\
                Arr100* (  rx)*(1.-ry)*(1.-rz)+\
                Arr101* (  rx)*(1.-ry)*(  rz)+\
                Arr110* (  rx)*(  ry)*(1.-rz)+\
                Arr111* (  rx)*(  ry)*(  rz);
    aget[0] = Aget;
    
    //return Aget;
}



__global__ void TraceBline(float *Bx,float *By,float *Bz,\
    float *inp_x,float *inp_y, float *inp_z,\
    float *out_x,float *out_y, float *out_z, unsigned long long N){

        unsigned long long x = blockIdx.x * blockDim.x + threadIdx.x;
        unsigned long long y = blockIdx.y * blockDim.y + threadIdx.y; 
        unsigned long long idx_cur = (gridDim.x*blockDim.x) * y + x;     
        if (idx_cur<N){
            // main procedure of B-line tracking 
        }
}
