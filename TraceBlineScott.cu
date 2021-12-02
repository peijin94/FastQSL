/**
*    @file TraceBline.cu
*    @author Peijin Zhang
*    The kernel of the Q-factor computation: tracing the magnetic field
*/
//#include <corecrt_math.h>
#include <math.h>
#include <stdio.h>
#include "helper_math.h"
#define M_PI 3.14159265    ///< Mathematical constant PI.
#define MAX_STEP_RATIO 2  ///< Maximum step length compared to box size.
#define INIT_DT 0.25 // default step length
#define N_STEP_LIM 10000 
#define TOL 1e-3 // toleranced error for each step [0.001~0.00001]

extern "C"{
#include "TraceBlineScott.cuh"
inline __device__ float lenVec3xyz(float xx,float yy,float zz){
    return sqrtf(xx*xx + yy*yy + zz*zz);}

inline __device__ float lenVec3(float3 a){
    return sqrtf(a.x*a.x + a.y*a.y + a.z*a.z);}

inline __device__ float dot3(float3 a, float3 b)
{  return a.x*b.x + a.y*b.y + a.z*b.z;}

inline __device__ double ddot3(float3 a, float3 b)
{  return (double)a.x*b.x + (double)a.y*b.y + (double)a.z*b.z;}

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

inline __device__ float3 Interp3dxyzn(float *Arr_x,float *Arr_y,float *Arr_z,\
    int3 AShapeN3, float3 inPoint_this,bool norm_flag){
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
    if(norm_flag){
        norm_B = sqrtf(res.x*res.x+res.y*res.y+res.z*res.z);
        res.x = res.x/norm_B;
        res.y = res.y/norm_B;
        res.z = res.z/norm_B;}
    return res;
}


inline __device__ float9 grad_unit_vec_B(float *Arr_x,float *Arr_y,float *Arr_z,\
    int3 AShapeN3, float3 inPoint_this){
        float3 vp1,vp2;
        float3 BN1,BN2;
        float9 res;

        vp1=inPoint_this; vp2=inPoint_this; vp1.x-=0.001; vp2.x+=0.001;
        BN1 = Interp3dxyzn(Arr_x,Arr_y,Arr_z,AShapeN3, vp1,true);
        BN2 = Interp3dxyzn(Arr_x,Arr_y,Arr_z,AShapeN3, vp2,true);
        res.x = (BN2-BN1)/0.002;

        vp1=inPoint_this; vp2=inPoint_this; vp1.y-=0.001; vp2.y+=0.001;
        BN1 = Interp3dxyzn(Arr_x,Arr_y,Arr_z,AShapeN3, vp1,true);
        BN2 = Interp3dxyzn(Arr_x,Arr_y,Arr_z,AShapeN3, vp2,true);
        res.y = (BN2-BN1)/0.002;

        vp1=inPoint_this; vp2=inPoint_this; vp1.z-=0.001; vp2.z+=0.001;
        BN1 = Interp3dxyzn(Arr_x,Arr_y,Arr_z,AShapeN3, vp1,true);
        BN2 = Interp3dxyzn(Arr_x,Arr_y,Arr_z,AShapeN3, vp2,true);
        res.z = (BN2-BN1)/0.002;

        return res;
}
    


inline __device__ float9 f_scott( float *Arr_x,float *Arr_y,float *Arr_z,\
    int3 AShapeN3, float9 vec_in){

        float3 BN1;
        float9 Jac_B,vec_out;
        BN1   = Interp3dxyzn(Arr_x,Arr_y,Arr_z,AShapeN3, vec_in.x ,true);
        Jac_B = grad_unit_vec_B(Arr_x,Arr_y,Arr_z,AShapeN3, vec_in.x);

        vec_out.x = BN1;

        vec_out.y.x = dot(vec_in.y, make_float3(Jac_B.x.x,Jac_B.y.x,Jac_B.z.x));
        vec_out.y.y = dot(vec_in.y, make_float3(Jac_B.x.y,Jac_B.y.y,Jac_B.z.y));
        vec_out.y.z = dot(vec_in.y, make_float3(Jac_B.x.z,Jac_B.y.z,Jac_B.z.z)); 

        vec_out.z.x = dot(vec_in.z, make_float3(Jac_B.x.x,Jac_B.y.x,Jac_B.z.x));
        vec_out.z.y = dot(vec_in.z, make_float3(Jac_B.x.y,Jac_B.y.y,Jac_B.z.y));
        vec_out.z.z = dot(vec_in.z, make_float3(Jac_B.x.z,Jac_B.y.z,Jac_B.z.z)); 

        return vec_out;
}



/**
* Adaptive step-size integral scheme
*/
inline __device__ float RKF45_Scott(float *Bx,float *By,float *Bz,\
    int3 BshapeN3, float9 vec_in, float9 *vec_out, float s_len){
    float9 k1,k2,k3,k4,k5,k6;
    float err_step; 

    // parameters of the Butcher tableau
    //float c2,c3,c4,c5,c6;
    // c need  to be included only when f' depends on x
    float a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65;
    float b1,b3,b4,b5,b6;
    float ce1,ce3,ce4,ce5,ce6;
    a21=1./4.;
    a31=3./32.;      a32=9./32.;
    a41=1932./2197.; a42=-7200./2197.; a43= 7296./2197.;
    a51=439./216.;   a52= -8.;          a53=  3680./513.;   a54=-845./4104.;
    a61= -8./27.;    a62= 2.;           a63= -3544./2565.;  a64= 1859./4104.; a65= -11./40.;
    b1 = 16./135.;    b3 = 6656./12825.;   b4 = 28561./56430.; b5 = -9./50.;  b6 = 2./55.; 
    //b2 =0.; 
    //bb1= 25./216.;   bb2=0.;  bb3 = 1408./2565.; bb4 = 2197./4104.;  bb5 = -1./5.;
    ce1 = 1./360.;     ce3 = -128./4275.;    ce4 = -2197./75240.;   ce5 = 1./50.;   ce6 = 2./55.;
    
    k1 = s_len*f_scott(Bx,By,Bz,BshapeN3,vec_in);
    k2 = s_len*f_scott(Bx,By,Bz,BshapeN3,vec_in+ (a21*k1));
    k3 = s_len*f_scott(Bx,By,Bz,BshapeN3,vec_in+ (a31*k1+ a32*k2));
    k4 = s_len*f_scott(Bx,By,Bz,BshapeN3,vec_in+ (a41*k1+ a42*k2+ a43*k3));
    k5 = s_len*f_scott(Bx,By,Bz,BshapeN3,vec_in+ (a51*k1+ a52*k2+ a53*k3+ a54*k4));
    k6 = s_len*f_scott(Bx,By,Bz,BshapeN3,vec_in+ (a61*k1+ a62*k2+ a63*k3+ a64*k4+ a65*k5));
    
    vec_out[0] = vec_in+ (b1*k1  +b3*k3+  b4*k4  +b5*k5  +b6*k6); //b2=0
    //P_b = P0+ (bb1*k1 +bb2*k2 +bb3*k3+ bb4*k4 +bb5*k5);

    //err_step = lenVec3(P_a-P_b);
    err_step = lenVec3(ce1*k1.x+ce3*k3.x+ce4*k4.x+ce5*k5.x+ce6*k6.x);
    return err_step;
}


inline __device__ float selectFloat3xyz(float3 a, int dim){
    float res;
    switch(dim) {
        case 0  : res = a.x;  break;
        case 1  : res = a.y;  break;
        case 2  : res = a.z;  break;}
    return res;
}
inline __device__ float selectInt3xyz(int3 a, int dim){
    int res;
    switch(dim) {
        case 0  : res = a.x;  break;
        case 1  : res = a.y;  break;
        case 2  : res = a.z;  break;}
    return res;
}

inline __device__ int checkFlag(int3 BshapeN3, float3 P_cur){
    // check current status
    int flag_res = 42; // 42 means un-categorized
    // flag=0 means inside running box
    if (P_cur.x>=0. &P_cur.y>0. &P_cur.z>=0. &  \
        P_cur.x<=BshapeN3.x-1. &P_cur.y<=BshapeN3.y-1. & P_cur.z<=BshapeN3.z-1. ){flag_res=0;} 
    else{ // ouside
        if (P_cur.x< 0.             &P_cur.y>=0.            &P_cur.z>=0. &  \
            P_cur.x< BshapeN3.x-1.  &P_cur.y< BshapeN3.y-1. &P_cur.z< BshapeN3.z-1.   ){flag_res=1;} // x min 
        if (P_cur.x>=0.             &P_cur.y>=0.            &P_cur.z>=0. &  \
            P_cur.x>=BshapeN3.x-1.  &P_cur.y< BshapeN3.y-1. &P_cur.z< BshapeN3.z-1.   ){flag_res=2;} // x max
        if (P_cur.x>=0.             &P_cur.y< 0.            &P_cur.z>=0. &  \
            P_cur.x< BshapeN3.x-1.  &P_cur.y< BshapeN3.y-1. &P_cur.z< BshapeN3.z-1.   ){flag_res=3;} // y min
        if (P_cur.x>=0.             &P_cur.y>=0.            &P_cur.z>=0. &  \
            P_cur.x< BshapeN3.x-1.  &P_cur.y>=BshapeN3.y-1. &P_cur.z< BshapeN3.z-1.   ){flag_res=4;} // y max 
        if (P_cur.x>=0.             &P_cur.y>=0.            &P_cur.z< 0. &  \
            P_cur.x< BshapeN3.x-1.  &P_cur.y< BshapeN3.y-1. &P_cur.z< BshapeN3.z-1.   ){flag_res=5;} // z min
        if (P_cur.x>=0.             &P_cur.y>=0.            &P_cur.z>=0. &  \
            P_cur.x< BshapeN3.x-1.  &P_cur.y< BshapeN3.y-1. &P_cur.z>=BshapeN3.z-1.   ){flag_res=6;} // z max
    }
    // dim = int((flag-1)/2)
    return flag_res;
}

__device__ float9 TraceBlineScott(float *Bx,float *By,float *Bz,int3 BshapeN3,\
    float *curB_x, float *curB_y,  float *curB_z,double *twist_this,bool *curB_flag,\
    float *P_0, float *ncross_dir, float s_len, int *flag_start, int *flag_end,\
    float *q_0, float *q_perp, double *len_this, float tol_coef){
        unsigned long step_count = 0;
        //unsigned long step_lim = (MAX_STEP_RATIO*(BshapeN3.x+BshapeN3.y+BshapeN3.z));
        double len_lim = (MAX_STEP_RATIO*1.0*(BshapeN3.x+BshapeN3.y+BshapeN3.z));
        float err_step;
        double len_record=0;
        double twist=0;
        float3 PP1; // point coordinate
        float3 B_P1,B_P2,B_Pstart, ncross_dir3,cur_P1;
        float3 Bv_s,Bv_e;
        float3 v0,vs,ve,u0,us,ue,vs1,ve1,us1,ue1;
        int dir_sign,dim_start,dim_end;
        float9 vec_tmp,vec_a,vec_b,vec_s,vec_e; // vec_a vec_b step forward
        // vec_s vec_e start end

        float scale,tol_this,dL,Nb0;
        float p_mid, p1,p2; // for linear interpolation
        int flag_this;
        int dim_out;
        float4 P_tmp;
        bool z0flag;

        float9 posi_vec;

        float BN_s,BN_e;


        PP1=make_float3(P_0[0],P_0[1],P_0[2]);
        ncross_dir3=make_float3(ncross_dir[0],ncross_dir[1],ncross_dir[2]);
        if(PP1.z<-1e-8){return posi_vec;} // quit if under 0

        B_Pstart = Interp3dxyzn(Bx,By,Bz,BshapeN3,PP1,true); // b0
        Nb0 = length(Interp3dxyzn(Bx,By,Bz,BshapeN3,PP1,false));
        // tol settings

        if(PP1.z<3){    
            if (fabsf(dot3(B_Pstart,ncross_dir3))<=0.05){tol_this=TOL/8e3;}
            else {tol_this=TOL*powf(fabsf(dot3(B_Pstart,ncross_dir3)),3);}
        }
        else{tol_this=TOL;}
        
        tol_this=tol_this*tol_coef;
        
        v0 = make_float3(B_Pstart.y,-B_Pstart.x,0.0);
        v0 = normalize( v0-dot(v0,B_Pstart)*B_Pstart );

        u0.x = B_Pstart.y * v0.z - B_Pstart.z * v0.y;
        u0.y = B_Pstart.z * v0.x - B_Pstart.x * v0.z;
        u0.z = B_Pstart.x * v0.y - B_Pstart.y * v0.x;
        u0 = normalize(u0) ;

        z0flag = (PP1.z)<1e-6;
        len_record=0;
        // forward(+1) and backward(-1) direction
        for(dir_sign=-1;dir_sign<2;dir_sign+=2){
            vec_a.x=PP1;
            vec_a.y=u0;
            vec_a.z=v0;

            if(z0flag){ // start from z0
                if(B_Pstart.z*dir_sign < 0){ // downward    
                    if(dir_sign==-1) {vec_s=vec_a; flag_start[0]=5;}
                    else             {vec_e=vec_a; flag_end[0]=5;}
                    continue;
                }
            }   
            
            dL = 0;
            s_len = abs(s_len); // keep positive
            flag_this = 0;  // start from flag=0
            while ( (flag_this==0) & (len_record<len_lim) & (step_count<N_STEP_LIM)){
                
                // forward one step and return error
                err_step = RKF45_Scott(Bx, By, Bz, BshapeN3, vec_a,  &vec_tmp, s_len*dir_sign);

                vec_b=vec_tmp;
                scale = powf(tol_this/err_step/11.09,0.2);
                //out test
                if (scale<0.618){s_len = s_len*0.618;// redo RK45 when the error is too large
                    continue; }
                s_len = s_len*scale;
                if (s_len>100.)  {s_len=100.;} // upper limit of the step size
                if (s_len<1./10.){s_len=1./10.;} //lower limit of the step size
                //len_record = len_record+s_len;
                dL = lenVec3(vec_b.x-vec_a.x);
                len_record = len_record+dL;
                if (curB_flag[0]){
                    cur_P1 = Interp3dxyzn(curB_x,curB_y,curB_z,BshapeN3,vec_b.x,false);
                    B_P1 = Interp3dxyzn(Bx,By,Bz,BshapeN3,vec_b.x,false);
                    twist = twist+dot3(cur_P1,B_P1)/dot3(B_P1,B_P1)/4.0/M_PI*dL;
                }
    
                flag_this = checkFlag(BshapeN3,vec_b.x);  // check status
                if (flag_this>0){ // out of box
                    len_record = len_record-dL; // reverse step len
                    if (curB_flag[0]){// reverse twist
                        twist = twist-dot3(cur_P1,B_P1)/dot3(B_P1,B_P1)/4.0/M_PI*dL;
                    }
                    if (flag_this<=6){ // step out from surface
                        
                        // linear estimation
                        dim_out = int((flag_this-1)/2);
                        p1 = selectFloat3xyz(vec_a.x,dim_out);
                        p2 = selectFloat3xyz(vec_b.x,dim_out);
    
                        if (fabsf(p1-p2)>0){
                            if (flag_this%2==1) {p_mid=0;} // step out from min surface
                            else                {p_mid=float(selectInt3xyz(BshapeN3,dim_out));} // step out from max surface
                            
                            if(dir_sign==-1) {vec_s=(vec_a* (p2-p_mid) + vec_b* (p_mid-p1))/(p2-p1); flag_start[0]=flag_this;}
                            else             {vec_e=(vec_a* (p2-p_mid) + vec_b* (p_mid-p1))/(p2-p1); flag_end[0]  =flag_this;}
                            len_record = len_record+fabsf(p_mid-p1)/(1e-4+fabsf(selectFloat3xyz(B_P1,dim_out)));    
                        }
                        else{
                            if(dir_sign==-1) {vec_s=vec_a; flag_start[0]=flag_this;}
                            else             {vec_e=vec_a; flag_end[0]=flag_this;  }
                        }
    
                        if (curB_flag[0]){
                            twist = twist+dot3(cur_P1,B_P1)/dot3(B_P1,B_P1)/4.0/M_PI \
                            *fabsf(p_mid-p1)/(1e-4+fabsf(selectFloat3xyz(B_P1,dim_out)));
                        }    
                    }
                    else{ // not important situation
                        if(dir_sign==-1) {vec_s=vec_a; flag_start[0]=flag_this;}
                        else             {vec_e=vec_a; flag_end[0]=flag_this;  }
                        q_0[0]=0;
                        q_perp[0]=0;
                        return posi_vec;
                    }
                }
                vec_a=vec_b;
                step_count=step_count+1;
            }
        }

        
        
        Bv_s = Interp3dxyzn(Bx,By,Bz,BshapeN3,vec_s.x,false);
        dim_start = int((flag_start[0]-1)/2);
        BN_s = selectFloat3xyz(Bv_s,dim_start);
        us = vec_s.y;
        vs = vec_s.z;
        
        Bv_e = Interp3dxyzn(Bx,By,Bz,BshapeN3,vec_e.x,false);
        dim_end = int((flag_end[0]-1)/2);
        BN_e = selectFloat3xyz(Bv_e,dim_end);
        ue = vec_e.y;
        ve = vec_e.z;
    
        //printf("%d, [%f,%f,%f] , [%f], [%f,%f,%f] %d \n",dim_start,Bv_s.x,Bv_s.y,Bv_s.z,BN_s, vec_s.x.x,vec_s.x.y,vec_s.x.z,flag_start[0]);
        
        us1=us-selectFloat3xyz(us,dim_start)/selectFloat3xyz(Bv_s,dim_start)*Bv_s;
        vs1=vs-selectFloat3xyz(vs,dim_start)/selectFloat3xyz(Bv_s,dim_start)*Bv_s;
        ue1=ue-selectFloat3xyz(ue,dim_end)/selectFloat3xyz(Bv_e,dim_end)*Bv_e;
        ve1=ve-selectFloat3xyz(ve,dim_end)/selectFloat3xyz(Bv_e,dim_end)*Bv_e;


        q_0[0] =abs( dot(ue1,ue1)*dot(vs1,vs1)       \
                +    dot(us1,us1)*dot(ve1,ve1)       \
               - 2.0*dot(ue1,ve1)*dot(us1,vs1))/    \
                (abs( (Nb0*Nb0) / (BN_s*BN_e)));

        ue1=ue-dot(ue,Bv_e)/length(Bv_e)*(Bv_e/length(Bv_e));
        ve1=ve-dot(ve,Bv_e)/length(Bv_e)*(Bv_e/length(Bv_e));
        us1=us-dot(us,Bv_s)/length(Bv_s)*(Bv_s/length(Bv_s));
        vs1=vs-dot(vs,Bv_s)/length(Bv_s)*(Bv_s/length(Bv_s));

        q_perp[0] =abs( dot(ue1,ue1)*dot(vs1,vs1)   \
                  +     dot(us1,us1)*dot(ve1,ve1)   \
                  - 2.0*dot(ue1,ve1)*dot(us1,vs1))/ \
            ( (Nb0*Nb0) / (length(Bv_s)*length(Bv_e)));

        //printf("[%d][%f]:%f  :%f  :%f\n",step_count,P1[0],P1[1],P1[2]);
        len_this[0] = len_record;
        twist_this[0] = twist;

        posi_vec.x = vec_s.x; // start
        posi_vec.y = PP1;
        posi_vec.z = vec_e.x; // end

        return posi_vec;
    }

__global__ void TraceAllBlineScott(float *Bx,float *By,float *Bz,int *BshapeN,\
    float *curB_x, float *curB_y,  float *curB_z,double *twist,bool *curB_flag,\
    float *inp_x,float *inp_y, float *inp_z, float *inp_cross_dir,\
    float *start_x,float *start_y, float *start_z, int *flag_start_arr,\
    float *end_x,  float *end_y,   float *end_z,   int *flag_end_arr,\
    float *B_this_x,float *B_this_y, float *B_this_z,\
    float *q_0_arr,float *q_perp_arr,\
    float *s_len,unsigned long long *N,double *LineLen,float *tol_coef){
        
        unsigned long long x = blockIdx.x * blockDim.x + threadIdx.x;
        unsigned long long y = blockIdx.y * blockDim.y + threadIdx.y; 
        unsigned long long idx_cur,dim_all,works_per_thread,Bline_ID,line_idx;
        int3 BshapeN3 = make_int3(BshapeN[0],BshapeN[1],BshapeN[2]);

        dim_all = (gridDim.x*blockDim.x*gridDim.y*blockDim.y); // upper lim 8192 
        idx_cur = (gridDim.x*blockDim.x) * y + x;                     
        works_per_thread = N[0]/dim_all+1;
        
        double *twist_this = new double[1];
        float *P_0 = new float[3];
        float *P_out = new float[3];
        int *flag_cur = new int[1];
        double len_this;
        float9 posi_vec;
        float q_0_this;
        float q_perp_this;

        int flag_start,flag_end;

        for (line_idx=0; line_idx<works_per_thread; line_idx++){
            //Bline_ID = works_per_thread*idx_cur+line_idx;
            Bline_ID = idx_cur+line_idx*dim_all;
            if (Bline_ID<N[0]){
                //printf("  %llu ",Bline_ID);
                // forward
                P_0[0] = inp_x[Bline_ID];
                P_0[1] = inp_y[Bline_ID];
                P_0[2] = inp_z[Bline_ID]; 
                twist_this[0]=0;

                posi_vec = TraceBlineScott(Bx,By,Bz,BshapeN3,curB_x,curB_y,curB_z,twist_this,curB_flag,\
                    P_0, inp_cross_dir,s_len[0],&flag_start,&flag_end,\
                    &q_0_this, &q_perp_this, &len_this, tol_coef[0]);
            
                start_x[Bline_ID] = posi_vec.x.x;
                start_y[Bline_ID] = posi_vec.x.y;
                start_z[Bline_ID] = posi_vec.x.z;

                end_x[Bline_ID] = posi_vec.z.x;
                end_y[Bline_ID] = posi_vec.z.y;
                end_z[Bline_ID] = posi_vec.z.z;
                flag_start_arr[Bline_ID] = flag_start;
                flag_end_arr[Bline_ID] = flag_end;
                LineLen[Bline_ID] = len_this;
                if (curB_flag[0]){twist[Bline_ID]=twist_this[0];}

                B_this_x[Bline_ID] = Interp3d(Bx,BshapeN3,P_0[0],P_0[1],P_0[2]);
                B_this_y[Bline_ID] = Interp3d(By,BshapeN3,P_0[0],P_0[1],P_0[2]);
                B_this_z[Bline_ID] = Interp3d(Bz,BshapeN3,P_0[0],P_0[1],P_0[2]);

                q_0_arr[Bline_ID] = q_0_this;
                q_perp_arr[Bline_ID] = q_perp_this;
                
                //if (fabsf(B_this_x[Bline_ID]*inp_cross_dir[0]+B_this_y[Bline_ID]*inp_cross_dir[1]+B_this_z[Bline_ID]*inp_cross_dir[2])*100.\
                //  <lenVec3xyz(B_this_x[Bline_ID],B_this_y[Bline_ID],B_this_z[Bline_ID])){
                //    B_flag[Bline_ID] = 1;}
                //else{B_flag[Bline_ID] = 0;}
                //printf("flag***:  %d  %d\n",flag_cur[0],flag_start[Bline_ID]);
            }
        }
        delete[] twist_this;
        delete[] P_0;
        delete[] P_out;
        delete[] flag_cur;
}
}