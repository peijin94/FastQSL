#include <math.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265
#endif

__device__ float lenVec3(float xx,float yy,float zz){return sqrt(xx*xx+yy*yy+zz*zz);}

__device__ float get_Idx3d(float *Arr,long *AShapeN,long xIdx,long yIdx,long zIdx){
    return Arr[xIdx* AShapeN[1]*AShapeN[2]  +  yIdx* AShapeN[2]  +  zIdx];}

__device__ float Interp3d(float *Arr,long *AShapeN, \
    float inPoint_x, float inPoint_y, float inPoint_z){

    //algorithm [https://core.ac.uk/download/pdf/44386053.pdf]
    float rx,ry,rz; // ratio of the point
    float Arr000,Arr001,Arr010,Arr011,Arr100,Arr101,Arr110,Arr111;
    float Aget;
    int x_Idx,y_Idx,z_Idx;
    
    // handle out of boundary problem by extending
    inPoint_x = (inPoint_x>0) ? inPoint_x : 0;
    inPoint_y = (inPoint_y>0) ? inPoint_y : 0;
    inPoint_z = (inPoint_z>0) ? inPoint_z : 0;
    inPoint_x = (inPoint_x<AShapeN[0]) ? inPoint_x : AShapeN[0];
    inPoint_y = (inPoint_y<AShapeN[1]) ? inPoint_y : AShapeN[1];
    inPoint_z = (inPoint_z<AShapeN[2]) ? inPoint_z : AShapeN[2];

    // ratio of the points to adjacent grid
    rx = inPoint_x-floorf(inPoint_x);
    ry = inPoint_y-floorf(inPoint_y);
    rz = inPoint_z-floorf(inPoint_z);

    // index of point in the down-side
    x_Idx = __float2int_rd(inPoint_x);
    y_Idx = __float2int_rd(inPoint_y);
    z_Idx = __float2int_rd(inPoint_z);

    // grid boundary
    Arr000 = get_Idx3d(Arr,AShapeN, x_Idx  ,y_Idx  ,z_Idx  );
    Arr001 = get_Idx3d(Arr,AShapeN, x_Idx  ,y_Idx  ,z_Idx+1);
    Arr010 = get_Idx3d(Arr,AShapeN, x_Idx  ,y_Idx+1,z_Idx  );
    Arr011 = get_Idx3d(Arr,AShapeN, x_Idx  ,y_Idx+1,z_Idx+1);
    Arr100 = get_Idx3d(Arr,AShapeN, x_Idx+1,y_Idx  ,z_Idx  );
    Arr101 = get_Idx3d(Arr,AShapeN, x_Idx+1,y_Idx  ,z_Idx+1);
    Arr110 = get_Idx3d(Arr,AShapeN, x_Idx+1,y_Idx+1,z_Idx  );
    Arr111 = get_Idx3d(Arr,AShapeN, x_Idx+1,y_Idx+1,z_Idx+1);

    Aget =      Arr000* (1.-rx)*(1.-ry)*(1.-rz)+\
                Arr001* (1.-rx)*(1.-ry)*(  rz)+\
                Arr010* (1.-rx)*(  ry)*(1.-rz)+\
                Arr011* (1.-rx)*(  ry)*(  rz)+\
                Arr100* (  rx)*(1.-ry)*(1.-rz)+\
                Arr101* (  rx)*(1.-ry)*(  rz)+\
                Arr110* (  rx)*(  ry)*(1.-rz)+\
                Arr111* (  rx)*(  ry)*(  rz);
    
    return Aget;
}

__device__ void stepForward(float *Bx,float *By,float *Bz,float *BshapeN,\
        float *P_start, float *P_end, float s_len, int *flag){
    float Bx_cur,By_cur,Bz_cur,B0;
    Bx_cur  = Interp3d(Bx,BshapeN,P_start[0],P_start[1],P_start[2]);
    By_cur  = Interp3d(Bx,BshapeN,P_start[0],P_start[1],P_start[2]);
    Bz_cur  = Interp3d(Bx,BshapeN,P_start[0],P_start[1],P_start[2]);
    B0 = lenVec3(Bx_cur,By_cur,Bz_cur)
    P_end[0] = P_start[0]+s*Bx_cur/B0
    P_end[1] = P_start[1]+s*By_cur/B0
    P_end[2] = P_start[2]+s*Bz_cur/B0
}


__global__ void test_Idx3d(float *Arr,long *AShapeN, long *getIdx,float *res){
    res[0] = get_Idx3d(Arr,AShapeN,getIdx[0],getIdx[1],getIdx[2]);
}

__global__ void test_Interp3d(float *Arr,long *AShapeN, float *inPoint,float *res){
    res[0] = Interp3d(Arr,AShapeN,inPoint[0],inPoint[1],inPoint[2]);
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
