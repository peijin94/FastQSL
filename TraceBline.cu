#include <math.h>
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.14159265
#endif
#define MAX_STEP_RATIO 5

__device__ float lenVec3(float xx,float yy,float zz){
    return sqrtf(xx*xx + yy*yy + zz*zz);}

__device__ float get_Idx3d(float *Arr,int *AShapeN,int Idx0,int Idx1,int Idx2){
    //return Arr[Idx0* AShapeN[1]*AShapeN[2]  +  Idx1* AShapeN[2]  +  Idx2];
    return Arr[Idx2* AShapeN[1]*AShapeN[0]  +  Idx1* AShapeN[0]  +  Idx0];
}
// test of the indexing of function
__global__ void  test_Idx3d(float *Arr,int *AShapeN, int *getIdx,float *res){
    printf("idx : %d    %d    %d   ", getIdx[0],getIdx[1],getIdx[2]);
    res[0] = get_Idx3d(Arr,AShapeN,getIdx[0],getIdx[1],getIdx[2]);
}

__device__ float Interp3d(float *Arr,int *AShapeN, \
    float inPoint_0, float inPoint_1, float inPoint_2){

    //algorithm [https://core.ac.uk/download/pdf/44386053.pdf]
    float rx,ry,rz; // ratio of the point
    float Arr000,Arr001,Arr010,Arr011,Arr100,Arr101,Arr110,Arr111;
    float Aget;
    int x_Idx,y_Idx,z_Idx;

    // handle out of boundary problem by extending
    inPoint_0 = (inPoint_0>0) ? inPoint_0 : 0.00001;
    inPoint_1 = (inPoint_1>0) ? inPoint_1 : 0.00001;
    inPoint_2 = (inPoint_2>0) ? inPoint_2 : 0.00001;
    inPoint_0 = (inPoint_0<AShapeN[0]-1) ? inPoint_0 : ((float)AShapeN[0]-1-0.00001);
    inPoint_1 = (inPoint_1<AShapeN[1]-1) ? inPoint_1 : ((float)AShapeN[1]-1-0.00001);
    inPoint_2 = (inPoint_2<AShapeN[2]-1) ? inPoint_2 : ((float)AShapeN[2]-1-0.00001);

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
    float *P_start, float *P_end,float s_len){
    float Bx_cur,By_cur,Bz_cur,B0;
    Bx_cur  = Interp3d(Bx,BshapeN,P_start[0],P_start[1],P_start[2]);
    By_cur  = Interp3d(By,BshapeN,P_start[0],P_start[1],P_start[2]);
    Bz_cur  = Interp3d(Bz,BshapeN,P_start[0],P_start[1],P_start[2]);
    B0 = lenVec3(Bx_cur,By_cur,Bz_cur);
    P_end[0] = P_start[0] + s_len*Bx_cur/B0;
    P_end[1] = P_start[1] + s_len*By_cur/B0;
    P_end[2] = P_start[2] + s_len*Bz_cur/B0;
    //printf("Inr:%f  :%f  :%f\n",s_len*Bx_cur/B0,s_len*By_cur/B0,s_len*Bz_cur/B0);
}

__device__ void RK4(float *Bx,float *By,float *Bz,int *BshapeN,\
    float *P_start, float *P_end, float s_len){
    
        float Bx_k1,By_k1,Bz_k1,B0_k1;
        float Bx_k2,By_k2,Bz_k2,B0_k2;
        float Bx_k3,By_k3,Bz_k3,B0_k3;
        float Bx_k4,By_k4,Bz_k4,B0_k4;

        Bx_k1  = Interp3d(Bx,BshapeN,P_start[0],P_start[1],P_start[2]);
        By_k1  = Interp3d(By,BshapeN,P_start[0],P_start[1],P_start[2]);
        Bz_k1  = Interp3d(Bz,BshapeN,P_start[0],P_start[1],P_start[2]);
        B0_k1  = lenVec3(Bx_k1,By_k1,Bz_k1);
        
        Bx_k2  = Interp3d(Bx,BshapeN,P_start[0]+s_len*Bx_k1/B0_k1/2,P_start[1]+s_len*By_k1/B0_k1/2,P_start[2]+s_len*Bz_k1/B0_k1/2);
        By_k2  = Interp3d(By,BshapeN,P_start[0]+s_len*Bx_k1/B0_k1/2,P_start[1]+s_len*By_k1/B0_k1/2,P_start[2]+s_len*Bz_k1/B0_k1/2);
        Bz_k2  = Interp3d(Bz,BshapeN,P_start[0]+s_len*Bx_k1/B0_k1/2,P_start[1]+s_len*By_k1/B0_k1/2,P_start[2]+s_len*Bz_k1/B0_k1/2);
        B0_k2  = lenVec3(Bx_k2,By_k2,Bz_k2);
        
        Bx_k3  = Interp3d(Bx,BshapeN,P_start[0]+s_len*Bx_k2/B0_k2/2,P_start[1]+s_len*By_k2/B0_k2/2,P_start[2]+s_len*Bz_k2/B0_k2/2);
        By_k3  = Interp3d(By,BshapeN,P_start[0]+s_len*Bx_k2/B0_k2/2,P_start[1]+s_len*By_k2/B0_k2/2,P_start[2]+s_len*Bz_k2/B0_k2/2);
        Bz_k3  = Interp3d(Bz,BshapeN,P_start[0]+s_len*Bx_k2/B0_k2/2,P_start[1]+s_len*By_k2/B0_k2/2,P_start[2]+s_len*Bz_k2/B0_k2/2);
        B0_k3  = lenVec3(Bx_k3,By_k3,Bz_k3);
        
        Bx_k4  = Interp3d(Bx,BshapeN,P_start[0]+s_len*Bx_k3/B0_k3,P_start[1]+s_len*By_k3/B0_k3,P_start[2]+s_len*Bz_k3/B0_k3);
        By_k4  = Interp3d(By,BshapeN,P_start[0]+s_len*Bx_k3/B0_k3,P_start[1]+s_len*By_k3/B0_k3,P_start[2]+s_len*Bz_k3/B0_k3);
        Bz_k4  = Interp3d(Bz,BshapeN,P_start[0]+s_len*Bx_k3/B0_k3,P_start[1]+s_len*By_k3/B0_k3,P_start[2]+s_len*Bz_k3/B0_k3);
        B0_k4  = lenVec3(Bx_k4,By_k4,Bz_k4);

        P_end[0] = P_start[0] + (1./6.)* s_len*(Bx_k1/B0_k1 + 2.0*Bx_k2/B0_k2 + 2.0*Bx_k3/B0_k3 + Bx_k4/B0_k4 );
        P_end[1] = P_start[1] + (1./6.)* s_len*(By_k1/B0_k1 + 2.0*By_k2/B0_k2 + 2.0*By_k3/B0_k3 + By_k4/B0_k4 );
        P_end[2] = P_start[2] + (1./6.)* s_len*(Bz_k1/B0_k1 + 2.0*Bz_k2/B0_k2 + 2.0*Bz_k3/B0_k3 + Bz_k4/B0_k4 );   
}

// merge in main 
__device__ int checkFlag(int *BshapeN, float *P_cur){
    // check current status
    int flag_res = 42; // 42 means un-categorized
    // flag=0 means inside running box
    if (P_cur[0]>=0. &P_cur[1]>0. &P_cur[2]>=0. &  \
        P_cur[0]<=BshapeN[0]-1. &P_cur[1]<=BshapeN[1]-1. & P_cur[2]<=BshapeN[2]-1. ){flag_res=0;} 
    // flag=1 means outside box below (normal end of simulation)
    //if (P_cur[0]>=0.           &P_cur[1]>=0.           &P_cur[2]<=0. &  \
    //   P_cur[0]<BshapeN[0]-1. &P_cur[1]<BshapeN[1]-1. & P_cur[2]<BshapeN[2]-1. ){flag_res=1;} 
    //printf("flag:%d\t",flag_res);
    return flag_res;
}

__device__ void TraceBline(float *Bx,float *By,float *Bz,int *BshapeN,\
    float *P_0, float *P_out, float s_len, int *flag, unsigned long long *step_len_this){
        
        
        unsigned long step_count = 0;
        float ratio_s;
        int flag_this;
        float Bz_tmp, direction;
        
        float *P_tmp = new float[3];
        float *P1 = new float[3];
        float *P2 = new float[3];
        
        flag_this = 0;  // start from flag=0
        
        P1[0] = P_0[0]; 
        P1[1] = P_0[1];
        P1[2] = P_0[2];

        Bz_tmp  = Interp3d(Bz,BshapeN,P1[0],P1[1],P1[2]);
        
        if (Bz_tmp>0) {direction=1.0;}
        else{direction=-1.0;}
        //printf("Direction:%f\n",direction);

        while ( (flag_this==0) & (step_count<(MAX_STEP_RATIO*4*BshapeN[0]))){
            // trace Bline step by step
            RK4(Bx,By,Bz,BshapeN,P1, P2, -s_len);
            //stepForward(Bx,By,Bz,BshapeN,P1, P2, s_len*direction);
            //printf("Cur[%d] : %f  :%f  :%f\n",step_count,P2[0],P2[1],P2[2]);
            
            flag_this = checkFlag(BshapeN,P2);  // check status
            if (flag_this>0){
                if (flag_this==1){
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
            step_count = step_count+1;        
            
            //printf("[%d][%f]:%f  :%f  :%f\n",step_count,direction,P1[0],P1[1],P1[2]);
        }
        //printf("[%d][%f]:%f  :%f  :%f\n",step_count,direction,P1[0],P1[1],P1[2]);
        step_len_this[0] = step_count;

        flag[0] = flag_this;
        delete[] P_tmp;
        delete[] P1;
        delete[] P2;
    }


__global__ void test_Interp3d(float *Arr,int *AShapeN, float *inPoint,float *res){
    res[0] = Interp3d(Arr,AShapeN,inPoint[0],inPoint[1],inPoint[2]);
}


__global__ void TraceAllBline(float *Bx,float *By,float *Bz,int *BshapeN,\
    float *inp_x,float *inp_y, float *inp_z,\
    float *out_x,float *out_y, float *out_z,       float *s_len,\
    int *flag_out,unsigned long long *N,unsigned long long *stepLineLen){
        
        unsigned long long x = blockIdx.x * blockDim.x + threadIdx.x;
        unsigned long long y = blockIdx.y * blockDim.y + threadIdx.y; 
        unsigned long long idx_cur = (gridDim.x*blockDim.x) * y + x;                     
        float *P_0 = new float[3];
        float *P_out = new float[3];
        int *flag_cur = new int[1];
        unsigned long long *step_len_this = new unsigned long long[1];

        if (idx_cur<N[0]){
            // main procedure of B-line tracking 
            P_0[0] = inp_x[idx_cur];
            P_0[1] = inp_y[idx_cur];
            P_0[2] = inp_z[idx_cur]; 

            TraceBline(Bx,By,Bz,BshapeN,P_0, P_out, s_len[0], flag_cur,step_len_this);
            out_x[idx_cur] = P_out[0];
            out_y[idx_cur] = P_out[1];
            out_z[idx_cur] = P_out[2];
            flag_out[idx_cur] = flag_cur[0];
            stepLineLen[idx_cur] = step_len_this[0];
            //printf("[%d], %f, %f, %f\n", flag_out[idx_cur] ,out_x[idx_cur],out_y[idx_cur],out_z[idx_cur] );
        }
        delete[] P_0;
        delete[] P_out;
        delete[] flag_cur;
}

__global__ void TestMem(int *flag_out){

}


__global__ void TraceBline_test(float *Bx,float *By,float *Bz,int *BshapeN,\
    float *P_0, float *P_out, float *s_len, int *flag, unsigned long long *step_len_this,\
    float *xx,float *yy, float *zz){ 
        unsigned long step_count = 0;
        float ratio_s;
        int flag_this;
        float Bz_tmp, direction;
        
        float *P_tmp = new float[3];
        float *P1 = new float[3];
        float *P2 = new float[3];
        
        flag_this = 0;  // start from flag=0
        
        P1[0] = P_0[0]; 
        P1[1] = P_0[1];
        P1[2] = P_0[2];

        Bz_tmp  = Interp3d(Bz,BshapeN,P1[0],P1[1],P1[2]);
        
        if (Bz_tmp>0) {direction=1.0;}
        else{direction=-1.0;}
        //printf("Direction:%f\n",direction);

        while ( (flag_this==0) & (step_count<10000-1)){
            // trace Bline step by step
            stepForward(Bx,By,Bz,BshapeN,P1, P2, s_len[0]*direction);
            //printf("Cur[%d] : %f  :%f  :%f\n",step_count,P2[0],P2[1],P2[2]);
            
            flag_this = checkFlag(BshapeN,P2);  // check status
            if (flag_this>0){
                if (flag_this==1){
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
            xx[step_count] = P2[0] ;
            yy[step_count] = P2[1] ;
            zz[step_count] = P2[2] ;
            P_tmp[0] = P1[0];   P_tmp[1] = P1[1];   P_tmp[2] = P1[2];
               P1[0] = P2[0];      P1[1] = P2[1];      P1[2] = P2[2];
            step_count = step_count+1;        
        }
        printf("[%d][%f]:%f  :%f  :%f\n",step_count,direction,P1[0],P1[1],P1[2]);
        step_len_this[0] = step_count;

        flag[0] = flag_this;
        delete[] P_tmp;
        delete[] P1;
        delete[] P2;
    }
