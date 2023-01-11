/**
*    @file TraceBline.cu
*    @author Peijin Zhang
*    The kernel of the Q-factor computation: tracing the magnetic field
*/
#include <math.h>
#include <stdio.h>
//#include "TraceBline.cuh"
#define M_PI 3.14159265   ///< Mathematical constant PI.
#define MAX_STEP_RATIO 4  ///< Maximum step length compared to box size.
extern "C"{
__device__ float lenVec3(float xx,float yy,float zz){
    return sqrtf(xx*xx + yy*yy + zz*zz);}

__forceinline__ __device__ float dot3(float3 a, float3 b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

__forceinline__ __device__ float3 divide3(float3 a, float b)
{
  return make_float3(a.x/b , a.y/b , a.z/b);
}


__forceinline__ __device__ float get_Idx3d(float *Arr,int3 AShapeN3,int Idx0,int Idx1,int Idx2){
    //return Arr[Idx0* AShapeN[1]*AShapeN[2]  +  Idx1* AShapeN[2]  +  Idx2];
    return Arr[Idx2* AShapeN3.y*AShapeN3.x  +  Idx1*AShapeN3.x  +  Idx0];
}

__global__ void  test_Idx3d(float *Arr,int *AShapeN, int *getIdx,float *res){
    printf("idx : %d    %d    %d   ", getIdx[0],getIdx[1],getIdx[2]);
    int3 AShapeN3 = make_int3(AShapeN[0],AShapeN[1],AShapeN[2]);
    res[0] = get_Idx3d(Arr,AShapeN3,getIdx[0],getIdx[1],getIdx[2]);
}

__device__ float Interp3d(float *Arr,int3 AShapeN3, \
    float inPoint_0, float inPoint_1, float inPoint_2){
    //algorithm [https://core.ac.uk/download/pdf/44386053.pdf]
    float rx,ry,rz; // ratio of the point
    float Arr000,Arr001,Arr010,Arr011,Arr100,Arr101,Arr110,Arr111;
    float Aget;
    int x_Idx,y_Idx,z_Idx;

    // handle out of boundary problem by extending
    inPoint_0 = (inPoint_0>0) ? inPoint_0 : 0.0001;
    inPoint_1 = (inPoint_1>0) ? inPoint_1 : 0.0001;
    inPoint_2 = (inPoint_2>0) ? inPoint_2 : 0.0001;
    inPoint_0 = (inPoint_0<AShapeN3.x-1) ? inPoint_0 : ((float)AShapeN3.x-1-0.0001);
    inPoint_1 = (inPoint_1<AShapeN3.y-1) ? inPoint_1 : ((float)AShapeN3.y-1-0.0001);
    inPoint_2 = (inPoint_2<AShapeN3.z-1) ? inPoint_2 : ((float)AShapeN3.z-1-0.0001);

    // ratio of the points to adjacent grid
    rx = inPoint_0-floorf(inPoint_0);
    ry = inPoint_1-floorf(inPoint_1);
    rz = inPoint_2-floorf(inPoint_2);

    // index of point in the down-side
    x_Idx = __float2int_rd(inPoint_0);
    y_Idx = __float2int_rd(inPoint_1);
    z_Idx = __float2int_rd(inPoint_2);
    
    //printf("x,y,z %f %f %f\n",inPoint_0,inPoint_1,inPoint_2);
    //printf("x,y,z %d %d %d\n",x_Idx,y_Idx,z_Idx);

    // grid boundary
    Arr000 = get_Idx3d(Arr,AShapeN3, x_Idx  ,y_Idx  ,z_Idx  );
    Arr001 = get_Idx3d(Arr,AShapeN3, x_Idx  ,y_Idx  ,z_Idx+1);
    Arr010 = get_Idx3d(Arr,AShapeN3, x_Idx  ,y_Idx+1,z_Idx  );
    Arr011 = get_Idx3d(Arr,AShapeN3, x_Idx  ,y_Idx+1,z_Idx+1);
    Arr100 = get_Idx3d(Arr,AShapeN3, x_Idx+1,y_Idx  ,z_Idx  );
    Arr101 = get_Idx3d(Arr,AShapeN3, x_Idx+1,y_Idx  ,z_Idx+1);
    Arr110 = get_Idx3d(Arr,AShapeN3, x_Idx+1,y_Idx+1,z_Idx  );
    Arr111 = get_Idx3d(Arr,AShapeN3, x_Idx+1,y_Idx+1,z_Idx+1);

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

inline __device__ float3 Interp3dxyz(float *Arr_x,float *Arr_y,float *Arr_z,int3 AShapeN3, \
    float inPoint_0, float inPoint_1, float inPoint_2){

    //algorithm [https://core.ac.uk/download/pdf/44386053.pdf]
    float rx,ry,rz; // ratio of the point
    float Arr000_x,Arr001_x,Arr010_x,Arr011_x,Arr100_x,Arr101_x,Arr110_x,Arr111_x;
    float Arr000_y,Arr001_y,Arr010_y,Arr011_y,Arr100_y,Arr101_y,Arr110_y,Arr111_y;
    float Arr000_z,Arr001_z,Arr010_z,Arr011_z,Arr100_z,Arr101_z,Arr110_z,Arr111_z;
    float w000,  w001,  w010,  w011,  w100,  w101,  w110,  w111;
    int x_Idx,y_Idx,z_Idx;
    float3 res;
    // ratio of the points to adjacent grid
    rx = inPoint_0-floorf(inPoint_0);
    ry = inPoint_1-floorf(inPoint_1);
    rz = inPoint_2-floorf(inPoint_2);
    // index of point in the down-side
    x_Idx = __float2int_rd(inPoint_0);
    y_Idx = __float2int_rd(inPoint_1);
    z_Idx = __float2int_rd(inPoint_2);
    
    if (inPoint_0<= 0.){
        rx=0; x_Idx=0;}
    if (inPoint_0>= AShapeN3.x -1){
        rx=1; x_Idx=AShapeN3.x -2;}
    if (inPoint_1<= 0.){
        ry=0; y_Idx=0;}
    if (inPoint_1>= AShapeN3.y -1){
        ry=1; y_Idx=AShapeN3.y -2;}
    if (inPoint_2<= 0.){
        rz=0.; z_Idx=0;}
    if (inPoint_2>= AShapeN3.z -1){
        rz=1.; z_Idx=AShapeN3.z -2;}
    
    // calculate the weight of the point and use for three times
    w000 = (1.-rx)*(1.-ry)*(1.-rz);
    w001 = (1.-rx)*(1.-ry)*(   rz);
    w010 = (1.-rx)*(   ry)*(1.-rz);
    w011 = (1.-rx)*(   ry)*(   rz);
    w100 = (   rx)*(1.-ry)*(1.-rz);
    w101 = (   rx)*(1.-ry)*(   rz);
    w110 = (   rx)*(   ry)*(1.-rz);
    w111 = (   rx)*(   ry)*(   rz);

    // grid boundary
    Arr000_x = get_Idx3d(Arr_x,AShapeN3, x_Idx  ,y_Idx  ,z_Idx  );
    Arr001_x = get_Idx3d(Arr_x,AShapeN3, x_Idx  ,y_Idx  ,z_Idx+1);
    Arr010_x = get_Idx3d(Arr_x,AShapeN3, x_Idx  ,y_Idx+1,z_Idx  );
    Arr011_x = get_Idx3d(Arr_x,AShapeN3, x_Idx  ,y_Idx+1,z_Idx+1);
    Arr100_x = get_Idx3d(Arr_x,AShapeN3, x_Idx+1,y_Idx  ,z_Idx  );
    Arr101_x = get_Idx3d(Arr_x,AShapeN3, x_Idx+1,y_Idx  ,z_Idx+1);
    Arr110_x = get_Idx3d(Arr_x,AShapeN3, x_Idx+1,y_Idx+1,z_Idx  );
    Arr111_x = get_Idx3d(Arr_x,AShapeN3, x_Idx+1,y_Idx+1,z_Idx+1);
    
    Arr000_y = get_Idx3d(Arr_y,AShapeN3, x_Idx  ,y_Idx  ,z_Idx  );
    Arr001_y = get_Idx3d(Arr_y,AShapeN3, x_Idx  ,y_Idx  ,z_Idx+1);
    Arr010_y = get_Idx3d(Arr_y,AShapeN3, x_Idx  ,y_Idx+1,z_Idx  );
    Arr011_y = get_Idx3d(Arr_y,AShapeN3, x_Idx  ,y_Idx+1,z_Idx+1);
    Arr100_y = get_Idx3d(Arr_y,AShapeN3, x_Idx+1,y_Idx  ,z_Idx  );
    Arr101_y = get_Idx3d(Arr_y,AShapeN3, x_Idx+1,y_Idx  ,z_Idx+1);
    Arr110_y = get_Idx3d(Arr_y,AShapeN3, x_Idx+1,y_Idx+1,z_Idx  );
    Arr111_y = get_Idx3d(Arr_y,AShapeN3, x_Idx+1,y_Idx+1,z_Idx+1);
    
    Arr000_z = get_Idx3d(Arr_z,AShapeN3, x_Idx  ,y_Idx  ,z_Idx  );
    Arr001_z = get_Idx3d(Arr_z,AShapeN3, x_Idx  ,y_Idx  ,z_Idx+1);
    Arr010_z = get_Idx3d(Arr_z,AShapeN3, x_Idx  ,y_Idx+1,z_Idx  );
    Arr011_z = get_Idx3d(Arr_z,AShapeN3, x_Idx  ,y_Idx+1,z_Idx+1);
    Arr100_z = get_Idx3d(Arr_z,AShapeN3, x_Idx+1,y_Idx  ,z_Idx  );
    Arr101_z = get_Idx3d(Arr_z,AShapeN3, x_Idx+1,y_Idx  ,z_Idx+1);
    Arr110_z = get_Idx3d(Arr_z,AShapeN3, x_Idx+1,y_Idx+1,z_Idx  );
    Arr111_z = get_Idx3d(Arr_z,AShapeN3, x_Idx+1,y_Idx+1,z_Idx+1);

    res.x=Arr000_x* w000+\
            Arr001_x* w001+\
            Arr010_x* w010+\
            Arr011_x* w011+\
            Arr100_x* w100+\
            Arr101_x* w101+\
            Arr110_x* w110+\
            Arr111_x* w111;
    res.y=Arr000_y* w000+\
            Arr001_y* w001+\
            Arr010_y* w010+\
            Arr011_y* w011+\
            Arr100_y* w100+\
            Arr101_y* w101+\
            Arr110_y* w110+\
            Arr111_y* w111;
    res.z=Arr000_z* w000+\
            Arr001_z* w001+\
            Arr010_z* w010+\
            Arr011_z* w011+\
            Arr100_z* w100+\
            Arr101_z* w101+\
            Arr110_z* w110+\
            Arr111_z* w111;

    return res;
}

__global__ void test_Interp3dxyz(float *Arr_x,float *Arr_y,float *Arr_z,int *AShapeN, \
    float *inPoint_0, float *inPoint_1, float *inPoint_2){
        float3 res;
        int3 shapeshape = make_int3(AShapeN[0],AShapeN[1],AShapeN[2]);
        res = Interp3dxyz(Arr_x,Arr_y,Arr_z,shapeshape,inPoint_0[0],inPoint_1[0],inPoint_2[0]);
    }

__device__ void stepForward(float *Bx,float *By,float *Bz,int3 BshapeN3,\
    float *P_start, float *P_end,float s_len){
    float Bx_cur,By_cur,Bz_cur,B0;
    Bx_cur  = Interp3d(Bx,BshapeN3,P_start[0],P_start[1],P_start[2]);
    By_cur  = Interp3d(By,BshapeN3,P_start[0],P_start[1],P_start[2]);
    Bz_cur  = Interp3d(Bz,BshapeN3,P_start[0],P_start[1],P_start[2]);
    B0 = lenVec3(Bx_cur,By_cur,Bz_cur);
    P_end[0] = P_start[0] + s_len*Bx_cur/B0;
    P_end[1] = P_start[1] + s_len*By_cur/B0;
    P_end[2] = P_start[2] + s_len*Bz_cur/B0;
    //printf("Inr:%f  :%f  :%f\n",s_len*Bx_cur/B0,s_len*By_cur/B0,s_len*Bz_cur/B0);
}

/**
* 4-stage Runge-Kutta integral function
*/
inline __device__ float3 RK4(float *Bx,float *By,float *Bz,int3 BshapeN3, float3 P0, float s_len){
    float B0_k1,B0_k2,B0_k3,B0_k4;
    float3 Bk1,Bk2,Bk3,Bk4,P_end;

    Bk1 = Interp3dxyz(Bx,By,Bz,BshapeN3,P0.x,P0.y,P0.z);
    B0_k1  = sqrtf(dot3(Bk1,Bk1));
    Bk2 = Interp3dxyz(Bx,By,Bz,BshapeN3,P0.x+s_len*Bk1.x/B0_k1/2.,P0.y+s_len*Bk1.y/B0_k1/2.,P0.z+s_len*Bk1.z/B0_k1/2.);
    B0_k2  = sqrtf(dot3(Bk2,Bk2));
    Bk3  = Interp3dxyz(Bx,By,Bz,BshapeN3,P0.x+s_len*Bk2.x/B0_k2/2.,P0.y+s_len*Bk2.y/B0_k2/2.,P0.z+s_len*Bk2.z/B0_k2/2.);
    B0_k3  = sqrtf(dot3(Bk3,Bk3));
    Bk4 = Interp3dxyz(Bx,By,Bz,BshapeN3,P0.x+s_len*Bk3.x/B0_k3,P0.y+s_len*Bk3.y/B0_k3,P0.z+s_len*Bk3.z/B0_k3);
    B0_k4  = sqrtf(dot3(Bk4,Bk4));

    P_end.x = P0.x + (1./6.)* s_len*( Bk1.x/B0_k1 + 2.0*Bk2.x/B0_k2 + 2.0*Bk3.x/B0_k3 + Bk4.x/B0_k4);
    P_end.y = P0.y + (1./6.)* s_len*( Bk1.y/B0_k1 + 2.0*Bk2.y/B0_k2 + 2.0*Bk3.y/B0_k3 + Bk4.y/B0_k4);
    P_end.z = P0.z + (1./6.)* s_len*( Bk1.z/B0_k1 + 2.0*Bk2.z/B0_k2 + 2.0*Bk3.z/B0_k3 + Bk4.z/B0_k4);   
    return P_end;
}


inline __device__ float3 RK4_boundary(float *Bx,float *By,float *Bz,int3 BshapeN3, float3 P0, float s_len,int b_dim){
    float B0_k1,B0_k2,B0_k3,B0_k4;
    float3 Bk1,Bk2,Bk3,Bk4,P_end,k1,k2,k3,k4;
    Bk1 = Interp3dxyz(Bx,By,Bz,BshapeN3,P0.x,P0.y,P0.z);
    
    switch(b_dim) {
        case 0  : B0_k1 = Bk1.x;  break;
        case 1  : B0_k1 = Bk1.y;  break;
        case 2  : B0_k1 = Bk1.z;  break;}
    k1=divide3(Bk1,B0_k1);
    k1.z=1.;
    Bk2 = Interp3dxyz(Bx,By,Bz,BshapeN3,P0.x+s_len*k1.x/2.,P0.y+s_len*k1.y/2.,P0.z+s_len*k1.z/2.);
    switch(b_dim) {
        case 0  : B0_k2 = Bk2.x;  break;
        case 1  : B0_k2 = Bk2.y;  break;
        case 2  : B0_k2 = Bk2.z;  break;}
    k2=divide3(Bk2,B0_k2);
    k2.z=1.;
    Bk3  = Interp3dxyz(Bx,By,Bz,BshapeN3,P0.x+s_len*k2.x/2.,P0.y+s_len*k2.y/2.,P0.z+s_len*k2.z/2.);
    switch(b_dim) {
        case 0  : B0_k3 = Bk3.x;  break;
        case 1  : B0_k3 = Bk3.y;  break;
        case 2  : B0_k3 = Bk3.z;  break;}
    k3=divide3(Bk3,B0_k3);
    k3.z=1.;
    Bk4 = Interp3dxyz(Bx,By,Bz,BshapeN3,P0.x+s_len*k3.x,P0.y+s_len*k3.y,P0.z+s_len*k3.z);
    switch(b_dim) {
        case 0  : B0_k4 = Bk4.x;  break;
        case 1  : B0_k4 = Bk4.y;  break;
        case 2  : B0_k4 = Bk4.z;  break;}    
    k4=divide3(Bk4,B0_k4);
    k4.z=1.;
    P_end.x = P0.x + (1./6.)* s_len*( k1.x + 2.0*k2.x + 2.0*k3.x + k4.x);
    P_end.y = P0.y + (1./6.)* s_len*( k1.y + 2.0*k2.y + 2.0*k3.y + k4.y);
    P_end.z = P0.z + (1./6.)* s_len*( k1.z + 2.0*k2.z + 2.0*k3.z + k4.z);   
    return P_end;
}

// merge in main 
inline __device__ int checkFlag(int3 BshapeN3, float3 P_cur){
    // check current status
    int flag_res = 42; // 42 means un-categorized
    // flag=0 means inside running box
    if (P_cur.x>=0. &P_cur.y>0. &P_cur.z>=0. &  \
        P_cur.x<=BshapeN3.x-1. &P_cur.y<=BshapeN3.y-1. & P_cur.z<=BshapeN3.z-1. ){flag_res=0;} 
    // flag=1 means outside box below (normal end of simulation)
    if (P_cur.x>=0.           &P_cur.y>=0.           &P_cur.z<0. &  \
       P_cur.x<BshapeN3.x-1. &P_cur.y<BshapeN3.y-1. & P_cur.z<BshapeN3.z-1. ){flag_res=1;} 
    //printf("flag:%d\t",flag_res);
    return flag_res;
}

__device__ void TraceBline(float *Bx,float *By,float *Bz,int3 BshapeN3,\
    float *P_0, float *P_out, float s_len, int *flag, unsigned long long *step_len_this){
        
        unsigned long step_count = 0;
        unsigned long step_lim = (MAX_STEP_RATIO*2*(BshapeN3.x+BshapeN3.y+BshapeN3.z));
        float ratio_z;
        int flag_this;
        float Bz_tmp, direction;
        float3 PP1,PP2;
        float3 B_P1,B_P2;

        flag_this = 0;  // start from flag=0
        PP1=make_float3(P_0[0],P_0[1],P_0[2]);
        Bz_tmp  = Interp3d(Bz,BshapeN3,PP1.x,PP1.y,PP1.z);
        
        if (Bz_tmp>0) {direction=1.0;}
        else{direction=-1.0;}
        //printf("Direction:%f\n",direction);

        while ( (flag_this==0) & (step_count<step_lim)){
            // trace Bline step by step
            PP2 = RK4(Bx,By,Bz,BshapeN3,PP1, s_len*direction);
            //stepForward(Bx,By,Bz,BshapeN3,P1, P2, s_len*direction);
            //printf("Cur[%d] : %f  :%f  :%f\n",step_count,P2[0],P2[1],P2[2]);
            flag_this = checkFlag(BshapeN3,PP2);  // check status
            if (flag_this>0){
                if (flag_this==1){
                    // linear estimation
                    B_P1 = Interp3dxyz(Bx,By,Bz,BshapeN3,PP1.x,PP1.y,PP1.z);
                    B_P2 = Interp3dxyz(Bx,By,Bz,BshapeN3,PP2.x,PP2.y,PP2.z);
                    if (1| fabsf(B_P1.z)*20.< sqrtf(dot3(B_P1,B_P1)) | fabsf(B_P2.z)*20.< sqrtf(dot3(B_P2,B_P2)) ){
                            P_out[0] = (PP1.x* (-PP2.z) + PP2.x* (PP1.z))/(PP1.z-PP2.z);
                            P_out[1] = (PP1.y* (-PP2.z) + PP2.y* (PP1.z))/(PP1.z-PP2.z);
                            P_out[2] = (PP1.z* (-PP2.z) + PP2.z* (PP1.z))/(PP1.z-PP2.z);
                            flag_this=43;
                    }
                    else{
                        // rk4 to the surface
                        PP2 = RK4_boundary(Bx,By,Bz,BshapeN3,PP1,-PP1.z,2);
                        P_out[0] = PP2.x;  P_out[1] = PP2.y;  P_out[2] = PP2.z;
                    }
                }
                else{ // ignore
                    P_out[0] = PP2.x;  P_out[1] = PP2.y;  P_out[2] = PP2.z;
                }
            }
            PP1=PP2;
            step_count = step_count+1;        
            //printf("[%d][%f]:%f  :%f  :%f\n",step_count,direction,P1[0],P1[1],P1[2]);
        }
        //printf("[%d][%f]:%f  :%f  :%f\n",step_count,direction,P1[0],P1[1],P1[2]);
        step_len_this[0] = step_count;
        flag[0] = flag_this;
    }


__global__ void test_Interp3d(float *Arr,int *AShapeN, float *inPoint,float *res){
    int3 AShapeN3 = make_int3(AShapeN[0],AShapeN[1],AShapeN[2]);
    res[0] = Interp3d(Arr,AShapeN3,inPoint[0],inPoint[1],inPoint[2]);
}


__global__ void TraceAllBline(float *Bx,float *By,float *Bz,int *BshapeN,\
    float *inp_x,float *inp_y, float *inp_z,\
    float *out_x,float *out_y, float *out_z, float *s_len,\
    float *B_inp_x,float *B_inp_y, float *B_inp_z,\
    float *B_out_x,float *B_out_y, float *B_out_z,\
    int *flag_out,unsigned long long *N,unsigned long long *stepLineLen){
        
        unsigned long long x = blockIdx.x * blockDim.x + threadIdx.x;
        unsigned long long y = blockIdx.y * blockDim.y + threadIdx.y; 
        unsigned long long idx_cur,dim_all,works_per_thread,Bline_ID,line_idx;
        int3 BshapeN3 = make_int3(BshapeN[0],BshapeN[1],BshapeN[2]);

        dim_all = (gridDim.x*blockDim.x*gridDim.y*blockDim.y); //8192
        idx_cur = (gridDim.x*blockDim.x) * y + x;                     
        works_per_thread = N[0]/dim_all+1;
        //printf("  %llu ",works_per_thread);
        
        float *P_0 = new float[3];
        float *P_out = new float[3];
        int *flag_cur = new int[1];
        unsigned long long *step_len_this = new unsigned long long[1];

        for (line_idx=0; line_idx<works_per_thread; line_idx++){
            //Bline_ID = works_per_thread*idx_cur+line_idx;
            Bline_ID = idx_cur+line_idx*dim_all;
            if (Bline_ID<N[0]){
                //printf("  %llu ",Bline_ID);
                // main procedure of B-line tracking 
                P_0[0] = inp_x[Bline_ID];
                P_0[1] = inp_y[Bline_ID];
                P_0[2] = inp_z[Bline_ID]; 
                // record input B
                B_inp_x[Bline_ID] = Interp3d(Bx,BshapeN3,P_0[0],P_0[1],P_0[2]);
                B_inp_y[Bline_ID] = Interp3d(By,BshapeN3,P_0[0],P_0[1],P_0[2]);
                B_inp_z[Bline_ID] = Interp3d(Bz,BshapeN3,P_0[0],P_0[1],P_0[2]);
                TraceBline(Bx,By,Bz,BshapeN3,P_0, P_out, s_len[0], flag_cur,step_len_this);
                out_x[Bline_ID] = P_out[0];
                out_y[Bline_ID] = P_out[1];
                out_z[Bline_ID] = P_out[2];
                flag_out[Bline_ID] = flag_cur[0];
                stepLineLen[Bline_ID] = step_len_this[0];
                //printf("[%d], %f, %f, %f\n", flag_out[idx_cur] ,out_x[idx_cur],out_y[idx_cur],out_z[idx_cur] );
                // record output B
                B_out_x[Bline_ID] = Interp3d(Bx,BshapeN3,P_out[0],P_out[1],P_out[2]);
                B_out_y[Bline_ID] = Interp3d(By,BshapeN3,P_out[0],P_out[1],P_out[2]);
                B_out_z[Bline_ID] = Interp3d(Bz,BshapeN3,P_out[0],P_out[1],P_out[2]);
                
            }
        }
        
        delete[] P_0;
        delete[] P_out;
        delete[] flag_cur;
}

__global__ void TestMem(int *flag_out){

}
}