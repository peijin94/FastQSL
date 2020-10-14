#include <math.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265
#endif
#define MAX_STEP_RATIO 32

__device__ float lenVec3(float xx,float yy,float zz){return sqrt(xx*xx+yy*yy+zz*zz);}

__device__ float get_Idx3d(float *Arr,int *AShapeN,int xIdx,int yIdx,int zIdx){
    return Arr[xIdx* AShapeN[1]*AShapeN[2]  +  yIdx* AShapeN[2]  +  zIdx];}

__device__ float Interp3d(float *Arr,int *AShapeN, \
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

__device__ void stepForward(float *Bx,float *By,float *Bz,int *BshapeN,\
        float *P_start, float *P_end, float s_len, int *flag){
    float Bx_cur,By_cur,Bz_cur,B0;
    Bx_cur  = Interp3d(Bx,BshapeN,P_start[0],P_start[1],P_start[2]);
    By_cur  = Interp3d(Bx,BshapeN,P_start[0],P_start[1],P_start[2]);
    Bz_cur  = Interp3d(Bx,BshapeN,P_start[0],P_start[1],P_start[2]);
    B0 = lenVec3(Bx_cur,By_cur,Bz_cur);
    P_end[0] = P_start[0]+s_len*Bx_cur/B0;
    P_end[1] = P_start[1]+s_len*By_cur/B0;
    P_end[2] = P_start[2]+s_len*Bz_cur/B0;
}

__device__ void RK4(float *Bx,float *By,float *Bz,float *BshapeN,\
    float *P_start, float *P_end, float s_len, int *flag){}

__device__ int checkFlag(float *BshapeN, float *P_cur){
    // check current status
    int flag_res = 42; // 42 means un-categorized
    // flag=0 means inside running box
    if (P_cur[0]>0 &P_cur[1]>0 &P_cur[2]>0 &  \
        P_cur[0]<BshapeN[0] &P_cur[1]<BshapeN[1]&P_cur[2]<BshapeN[2] ){flag_res=0;} 
    // flag=1 means outside box below (normal end of simulation)
    if (P_cur[0]>0 &P_cur[1]>0 &P_cur[2]<0 &  \
        P_cur[0]<BshapeN[0] &P_cur[1]<BshapeN[1]&P_cur[2]<BshapeN[2] ){flag_res=1;} 

    return flag_res;
}

__device__ void TraceBline(float *Bx,float *By,float *Bz,int *BshapeN,\
    float *P_0, float *P_out, float s_len, int *flag){
        
        int step_count = 0;
        float ratio_s;
        float *P_tmp = new float[3];
        float *P1 = new float[3];
        float *P2 = new float[3];
        int *flag_this = new int[1];
        
        flag_this[0] = flag[0];  // start from flag=0
        
        P1[0] = P_0[0]; P1[1] = P_0[1]; P1[2] = P_0[2];
        
        while ( (flag_this[0]==0) & (step_count<(MAX_STEP_RATIO*BshapeN[0]))){
            // trace Bline step by step
            stepForward(Bx,By,Bz,BshapeN,P1, P2, s_len, flag_this);
            
            flag_this[0] = checkFlag(BshapeN,P2);  // check status
            if (flag_this[0]>0){
                if (flag_this[0]==1){
                    // linear estimation
                    ratio_s = P1[2]/(P1[2]-P2[2]);
                    P_out[0] = P1[0]* (1-ratio_s) + P1[0]* (ratio_s);
                    P_out[1] = P1[1]* (1-ratio_s) + P1[1]* (ratio_s);
                    P_out[2] = 0;
                }
                else{ // ignore
                    P_out[0] = P2[0];  P_out[1] = P2[1];  P_out[2] = P2[2];
                }
            }
            P_tmp[0] = P1[0];   P_tmp[1] = P1[1];   P_tmp[2] = P1[2];
               P1[0] = P2[0];      P1[1] = P2[1];      P1[2] = P2[2];
        }
        flag[0] = flag_this[0];
    }

__global__ void test_Idx3d(float *Arr,int *AShapeN, int *getIdx,float *res){
    res[0] = get_Idx3d(Arr,AShapeN,getIdx[0],getIdx[1],getIdx[2]);
}

__global__ void test_Interp3d(float *Arr,int *AShapeN, float *inPoint,float *res){
    res[0] = Interp3d(Arr,AShapeN,inPoint[0],inPoint[1],inPoint[2]);
}



__global__ void TraceAllBline(float *Bx,float *By,float *Bz,float *BshapeN,\
    float *inp_x,float *inp_y, float *inp_z,\
    float *out_x,float *out_y, float *out_z,\
    float s_len,\
    int *flag_out,unsigned long long N){
        
        unsigned long long x = blockIdx.x * blockDim.x + threadIdx.x;
        unsigned long long y = blockIdx.y * blockDim.y + threadIdx.y; 
        unsigned long long idx_cur = (gridDim.x*blockDim.x) * y + x;     
        
        float *P_0 = new float[3];
        float *P_out = new float[3];
        int *flag_cur = new int[1];

        P_0[0] = inp_x[idx_cur];
        P_0[1] = inp_y[idx_cur];
        P_0[2] = inp_z[idx_cur];

        if (idx_cur<N){
            x_cur = inp_x[idx_cur];
            y_cur = inp_y[idx_cur];
            z_cur = inp_z[idx_cur];
            // main procedure of B-line tracking 
            TraceBline(Bx,By,Bz,BshapeN,P_0, P_out, s_len, flag_cur);

            out_x[idx_cur]=P_out[0]
            out_y[idx_cur]=P_out[1]
            out_z[idx_cur]=P_out[2]
        }
}
