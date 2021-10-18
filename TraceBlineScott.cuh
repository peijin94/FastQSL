/**
*    @file TraceBlineAdaptive.cuh
*    @author Peijin Zhang
*    The kernel of the Q-factor computation: tracing the magnetic field
*/


/**
 * Obtain value from a 3D array, for a given point P.
 * @param[in] Arr  The pointer of the input array
 * @param[in] AShapeN3  The shape of the input array, in the form of int3{x,y,z}
 * @param[in] Idx0  The index of P in the first demention
 * @param[in] Idx1  The index of P in the second demention
 * @param[in] Idx2  The index of P in the third demention
 * @return The value of P in the array.
 */
 __forceinline__ __device__ float get_Idx3d(float *Arr,int3 AShapeN3,int Idx0,int Idx1,int Idx2);

 /**
* A test terminal of the indexing of function "get_Idx3d"
*/
__global__ void  test_Idx3d(float *Arr,int *AShapeN, int *getIdx,float *res);

/**
* Interpolation of a 3D array for a given point.
* Assuming that the coordinate is the same as the grid index of xyz
* The point out of the range is evaluated by fixed value extrapolation
* Ref: [https://core.ac.uk/download/pdf/44386053.pdf]
* @param[in] Arr The pointer to the 3D array for interpolation
* @param[in] AShapeN3 The shape of the array
* @param[in] inPoint_0 Position in x axis
* @param[in] inPoint_0 Position in y axis
* @param[in] inPoint_0 Position in z axis
* @return The interpolated value
*/
__device__ float Interp3d(float *Arr,int3 AShapeN3, \
    float inPoint_0, float inPoint_1, float inPoint_2);

/**
* Interpolation of normalized Bx,By,Bz at one time, the position for interpolation is the same for all three arrays.
* The result is anormalized according to B0
* @param[in] Arr_x Bx, array for interpolation
* @param[in] Arr_y By, array for interpolation
* @param[in] Arr_z Bz, array for interpolation
* @param[in] AShapeN3 The shape of the arrays,(should be the same)
* @param[in] inPoint_0 Position in x axis
* @param[in] inPoint_0 Position in y axis
* @param[in] inPoint_0 Position in z axis
* @return The interpolated value of Bx,By,Bz,w
*/
inline __device__ float3 Interp3dxyzn(float *Arr_x,float *Arr_y,float *Arr_z,int3 AShapeN3, float3 inPoint_this,bool norm_flag);

inline __device__ float3 RK4(float *Bx,float *By,float *Bz,int3 BshapeN3, float3 P0, float s_len);

inline __device__ float4 RKF45(float *Bx,float *By,float *Bz,int3 BshapeN3, float3 P0, float s_len);

inline __device__ float3 RK4_boundary(float *Bx,float *By,float *Bz,int3 BshapeN3, float3 P0, float s_len,int b_dim);

inline __device__ int checkFlag(int3 BshapeN3, float3 P_cur);

__device__ void TraceBlineAdap(float *Bx,float *By,float *Bz,int3 BshapeN3,\
    float *curB_x, float *curB_y,  float *curB_z,double *twist,bool *curB_flag,\
    float *P_0, float *P_out, float *ncross_dir, float s_len, int *flag, double *len_this,\
    float direction,float tol_coef);

__global__ void test_Interp3d(float *Arr,int *AShapeN, float *inPoint,float *res);


__global__ void TraceAllBline(float *Bx,float *By,float *Bz,int *BshapeN,\
    float *curB_x, float *curB_y,  float *curB_z,double *twist,bool *curB_flag,\
    float *inp_x,float *inp_y, float *inp_z, float *inp_cross_dir,\
    float *start_x,float *start_y, float *start_z, int *flag_start,\
    float *end_x,  float *end_y,   float *end_z,   int *flag_end,\
    float *B_this_x,float *B_this_y, float *B_this_z, int *B_flag,\
    float *B_start_x,float *B_start_y, float *B_start_z,\
    float *B_end_x,float *B_end_y, float *B_end_z,\
    float *s_len,unsigned long long *N,double *LineLen,float *tol_coef);

__global__ void TestMem(int *flag_out);