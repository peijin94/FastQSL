/**
*    @file TraceBline.cu
*    @author Peijin Zhang
*    The kernel of the Q-factor computation: tracing the magnetic field
*/
#include <math.h>
#include <stdio.h>
#include "helper_math.h"
#define M_PI 3.14159265   ///< Mathematical constant PI.
#define MAX_STEP_RATIO 4  ///< Maximum step length compared to box size.
extern "C"{
#include "TraceBlineAdaptive.cuh"
inline __device__ float lenVec3xyz(float xx,float yy,float zz){
    return sqrtf(xx*xx + yy*yy + zz*zz);}

inline __device__ float lenVec3(float3 a){
    return sqrtf(a.x*a.x + a.y*a.y + a.z*a.z);}

inline __device__ float dot3(float3 a, float3 b)
{  return a.x*b.x + a.y*b.y + a.z*b.z;}

inline __device__ float3 divide3(float3 a, float b)
{  return make_float3(a.x/b , a.y/b , a.z/b);}

inline __device__ float3 normalize3(float4 a)
{  return make_float3(a.x/a.w,  a.y/a.w,  a.z/a.w);}

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

inline __device__ float3 Interp3dxyzn(float *Arr_x,float *Arr_y,float *Arr_z,int3 AShapeN3, float3 inPoint_this){
    // normalized B interpolation
    //algorithm [https://core.ac.uk/download/pdf/44386053.pdf]
    float rx,ry,rz; // ratio of the point
    float Arr000_x,Arr001_x,Arr010_x,Arr011_x,Arr100_x,Arr101_x,Arr110_x,Arr111_x;
    float Arr000_y,Arr001_y,Arr010_y,Arr011_y,Arr100_y,Arr101_y,Arr110_y,Arr111_y;
    float Arr000_z,Arr001_z,Arr010_z,Arr011_z,Arr100_z,Arr101_z,Arr110_z,Arr111_z;
    float w000,  w001,  w010,  w011,  w100,  w101,  w110,  w111;
    int x_Idx,y_Idx,z_Idx;
    float norm_B;
    float3 res;

    float inPoint_0=inPoint_this.x;
    float inPoint_1=inPoint_this.y;
    float inPoint_2=inPoint_this.z;
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
    norm_B = sqrtf(res.x*res.x+res.y*res.y+res.z*res.z);
    res.x = res.x/norm_B;
    res.y = res.y/norm_B;
    res.z = res.z/norm_B;
    return res;
}


__global__ void test_Interp3dxyz(float *Arr_x,float *Arr_y,float *Arr_z,int *AShapeN, \
    float *inPoint_0, float *inPoint_1, float *inPoint_2){
        float3 res,pointthis;
        pointthis=make_float3(inPoint_0[0],inPoint_1[0],inPoint_2[0]);
        int3 shapeshape = make_int3(AShapeN[0],AShapeN[1],AShapeN[2]);
        res = Interp3dxyzn(Arr_x,Arr_y,Arr_z,shapeshape,pointthis);
    }

/**
* Adaptive step-size integral scheme
*/
inline __device__ float4 RKF45(float *Bx,float *By,float *Bz,int3 BshapeN3, float3 P0, float s_len){
    float3 k1,k2,k3,k4,k5,k6;
    float3 P_a,P_b;
    float4 res_end;
    float err_step; 

    // parameters of the Butcher tableau
    //float c2,c3,c4,c5,c6;
    // c need  to be included only when f' depends on x
    float a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65;
    float b1,b2,b3,b4,b5,b6;
    float bb1,bb2,bb3,bb4,bb5; 
    float ce1,ce3,ce4,ce5,ce6;
    //c2=1./4.;   c3 = 3./8.;  c4=12./13.; c5=1.;  c6=1./2.; 
    a21=1./4.;
    a31=3./32.;      a32=9./32.;
    a41=1932./2197.; a42=-7200./2197.; a43= 7296./2197.;
    a51=439./216.;   a52= -8.;          a53=  3680./513.;   a54=-845./4104.;
    a61= -8./27.;    a62= 2.;           a63= -3544./2565.;  a64= 1859./4104.; a65= -11./40.;
    b1 = 16./135.;   b2 =0.;  b3 = 6656./12825.;   b4 = 28561./56430.; b5 = -9./50.;  b6 = 2./55.; 
    bb1= 25./216.;   bb2=0.;  bb3 = 1408./2565.; bb4 = 2197./4104.;  bb5 = -1./5.;
    ce1 = 1./360.;     ce3 = -128./4275.;    ce4 = -2197./75240.;   ce5 = 1./50.;   ce6 = 2./55.;
    k1 = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0);
    k2 = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0+ s_len*(a21*k1));
    k3 = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0+ s_len*(a31*k1+ a32*k2));
    k4 = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0+ s_len*(a41*k1+ a42*k2+ a43*k3));
    k5 = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0+ s_len*(a51*k1+ a52*k2+ a53*k3+ a54*k4));
    k6 = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0+ s_len*(a61*k1+ a62*k2+ a63*k3+ a64*k4+ a65*k5));
    
    P_a = P0+ s_len*(b1*k1  +b2*k2  +b3*k3+  b4*k4  +b5*k5  +b6*k6);
    P_b = P0+ s_len*(bb1*k1 +bb2*k2 +bb3*k3+ bb4*k4 +bb5*k5);

    err_step = lenVec3(ce1*k1+ce3*k3+ce4*k4+ce5*k5+ce6*k6);
    res_end.x = P_a.x;
    res_end.y = P_a.y;
    res_end.z = P_a.z;
    res_end.w = err_step;
    return res_end;
}

/**
* 4-stage Runge-Kutta integral function
*/
inline __device__ float3 RK4(float *Bx,float *By,float *Bz,int3 BshapeN3, float3 P0, float s_len){
    float3 k1,k2,k3,k4,P_end;
    k1 = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0);
    k2 = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0+s_len*k1/2.);
    k3 = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0+s_len*k2/2.);
    k4 = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0+s_len*k3);

    P_end = P0+1./6.*s_len*(k1 +2.*k2 +2.*k3 +k4);
    return P_end;
}

inline __device__ float3 RK4_boundary(float *Bx,float *By,float *Bz,int3 BshapeN3, float3 P0, float s_len,int b_dim){
    float B0_k1,B0_k2,B0_k3,B0_k4;
    float3 Bk1,Bk2,Bk3,Bk4,P_end,k1,k2,k3,k4;
    Bk1 = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0);
    switch(b_dim) {
        case 0  : B0_k1 = Bk1.x;  break;
        case 1  : B0_k1 = Bk1.y;  break;
        case 2  : B0_k1 = Bk1.z;  break;}
    k1 = Bk1/B0_k1;
    k1.z=1.;
    Bk2 = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0+s_len*k1/2.);
    switch(b_dim) {
        case 0  : B0_k2 = Bk2.x;  break;
        case 1  : B0_k2 = Bk2.y;  break;
        case 2  : B0_k2 = Bk2.z;  break;}
    k2 = Bk2/B0_k2;
    Bk2.z=1.;
    Bk3  = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0+s_len*k2/2.);
    switch(b_dim) {
        case 0  : B0_k3 = Bk3.x;  break;
        case 1  : B0_k3 = Bk3.y;  break;
        case 2  : B0_k3 = Bk3.z;  break;}
    k3 = Bk3/B0_k3;
    Bk3.z=1.;
    Bk4 = Interp3dxyzn(Bx,By,Bz,BshapeN3,P0+s_len*k3);
    switch(b_dim) {
        case 0  : B0_k4 = Bk4.x;  break;
        case 1  : B0_k4 = Bk4.y;  break;
        case 2  : B0_k4 = Bk4.z;  break;}    
    k4 = Bk4/B0_k4;
    Bk4.z=1.;
    P_end = P0 + (1./6.)* s_len*( k1 + 2.0*k2 + 2.0*k3 + k4);
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
        float4 P_tmp;

        flag_this = 0;  // start from flag=0
        PP1=make_float3(P_0[0],P_0[1],P_0[2]);
        Bz_tmp  = Interp3d(Bz,BshapeN3,PP1.x,PP1.y,PP1.z);
        
        if (Bz_tmp>0) {direction=1.0;}
        else{direction=-1.0;}
        //printf("Direction:%f\n",direction);

        while ( (flag_this==0) & (step_count<step_lim)){
            // trace Bline step by step
            
            PP2 = RK4(Bx,By,Bz,BshapeN3,PP1, s_len*direction);
            
            // P_tmp = RKF45(Bx,By,Bz,BshapeN3,PP1, s_len*direction);
            // PP2 = make_float3(P_tmp.x,P_tmp.y,P_tmp.z);
            
            //stepForward(Bx,By,Bz,BshapeN3,P1, P2, s_len*direction);
            //printf("Cur[%d] : %f  :%f  :%f\n",step_count,P2[0],P2[1],P2[2]);
            
            flag_this = checkFlag(BshapeN3,PP2);  // check status
            if (flag_this>0){
                if (flag_this==1){
                    // linear estimation
                    B_P1 = Interp3dxyzn(Bx,By,Bz,BshapeN3,PP1);
                    B_P2 = Interp3dxyzn(Bx,By,Bz,BshapeN3,PP2);
                    if ( fabsf(B_P1.z)*50.< sqrtf(dot3(B_P1,B_P1)) | fabsf(B_P2.z)*50.< sqrtf(dot3(B_P2,B_P2)) ){
                            P_out[0] = (PP1.x* (-PP2.z) + PP2.x* (PP1.z))/(PP1.z-PP2.z);
                            P_out[1] = (PP1.y* (-PP2.z) + PP2.y* (PP1.z))/(PP1.z-PP2.z);
                            P_out[2] = (PP1.z* (-PP2.z) + PP2.z* (PP1.z))/(PP1.z-PP2.z);
                            flag_this=10;
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
        
        }
        //printf("[%d][%f]:%f  :%f  :%f\n",step_count,direction,P1[0],P1[1],P1[2]);
        step_len_this[0] = step_count;
        flag[0] = flag_this;
}
__global__ void testTraceBlineAdap(float *Bx,float *By,float *Bz,int* BshapeN3,\
    float *P_0, float *P_out, float* s_len, int *flag, unsigned long long *step_len_this){
        TraceBlineAdap(Bx,By,Bz, make_int3(BshapeN3[0],BshapeN3[1],BshapeN3[2]),\
         P_0, P_out, s_len[0], flag, step_len_this);
    }
__device__ void TraceBlineAdap(float *Bx,float *By,float *Bz,int3 BshapeN3,\
    float *P_0, float *P_out, float s_len, int *flag, unsigned long long *step_len_this){
        
        unsigned long step_count = 0;
        unsigned long step_lim = (MAX_STEP_RATIO*2*(BshapeN3.x+BshapeN3.y+BshapeN3.z));
        float ratio_z,scale;
        int flag_this;
        float Bz_tmp, direction;
        float3 PP1,PP2;
        float3 B_P1,B_P2;
        float4 P_tmp;
        
        float tol=0.000001; // position error each step

        flag_this = 0;  // start from flag=0
        PP1=make_float3(P_0[0],P_0[1],P_0[2]);
        Bz_tmp  = Interp3d(Bz,BshapeN3,PP1.x,PP1.y,PP1.z);
        
        if (Bz_tmp>0) {direction=1.0;}
        else{direction=-1.0;}
        //printf("Direction:%f\n",direction);

        while ( (flag_this==0) & (step_count<step_lim)){
            // trace Bline step by step
            
            P_tmp = RKF45(Bx,By,Bz,BshapeN3,PP1, s_len*direction);
            PP2 = make_float3(P_tmp.x,P_tmp.y,P_tmp.z);
            scale = 0.85*powf(tol/P_tmp.w,0.25);
            if (scale<0.5){ // redo RK45 when the error is too large
                s_len = s_len*0.5;
                continue;
            }
            s_len = s_len*scale;
            if (s_len>128.){s_len=128.;}
            if (s_len<1./8.){s_len=1./8.;}

            //printf("[%d]  %f  %f  %f  %f  %f\n",step_count,PP2.x,PP2.y,PP2.z,scale,s_len);
            
            flag_this = checkFlag(BshapeN3,PP2);  // check status
            if (flag_this>0){
                if (flag_this==1){
                    // linear estimation
                    B_P1 = Interp3dxyzn(Bx,By,Bz,BshapeN3,PP1);
                    B_P2 = Interp3dxyzn(Bx,By,Bz,BshapeN3,PP2);
                    if ( fabsf(B_P1.z)*50.< sqrtf(dot3(B_P1,B_P1)) | fabsf(B_P2.z)*50.< sqrtf(dot3(B_P2,B_P2)) ){
                            P_out[0] = (PP1.x* (-PP2.z) + PP2.x* (PP1.z))/(PP1.z-PP2.z);
                            P_out[1] = (PP1.y* (-PP2.z) + PP2.y* (PP1.z))/(PP1.z-PP2.z);
                            P_out[2] = (PP1.z* (-PP2.z) + PP2.z* (PP1.z))/(PP1.z-PP2.z);
                            flag_this=8;
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
                TraceBlineAdap(Bx,By,Bz,BshapeN3,P_0, P_out, s_len[0], flag_cur,step_len_this);
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