__device__ float lenVec3(float xx,float yy,float zz);

__forceinline__ __device__ float dot3(float3 a, float3 b);

__forceinline__ __device__ float3 divide3(float3 a, float b);
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
* Interpolation of Bx,By,Bz at one time, the position for interpolation is the same for all three arrays.
* @param[in] Arr_x Bx, array for interpolation
* @param[in] Arr_y By, array for interpolation
* @param[in] Arr_z Bz, array for interpolation
* @param[in] AShapeN3 The shape of the arrays,(should be the same)
* @param[in] inPoint_0 Position in x axis
* @param[in] inPoint_0 Position in y axis
* @param[in] inPoint_0 Position in z axis
* @return The interpolated value of Bx,By,Bz
*/
inline __device__ float3 Interp3dxyz(float *Arr_x,float *Arr_y,float *Arr_z,int3 AShapeN3, \
    float inPoint_0, float inPoint_1, float inPoint_2);

__device__ void stepForward(float *Bx,float *By,float *Bz,int3 BshapeN3,\
    float *P_start, float *P_end,float s_len);

inline __device__ float3 RK4(float *Bx,float *By,float *Bz,int3 BshapeN3, float3 P0, float s_len);

inline __device__ float3 RK4_boundary(float *Bx,float *By,float *Bz,int3 BshapeN3, float3 P0, float s_len,int b_dim);

inline __device__ int checkFlag(int3 BshapeN3, float3 P_cur);

__device__ void TraceBline(float *Bx,float *By,float *Bz,int3 BshapeN3,\
    float *P_0, float *P_out, float s_len, int *flag, unsigned long long *step_len_this);

__global__ void test_Interp3d(float *Arr,int *AShapeN, float *inPoint,float *res);

__global__ void TraceAllBline(float *Bx,float *By,float *Bz,int *BshapeN,\
    float *inp_x,float *inp_y, float *inp_z,\
    float *out_x,float *out_y, float *out_z, float *s_len,\
    float *B_inp_x,float *B_inp_y, float *B_inp_z,\
    float *B_out_x,float *B_out_y, float *B_out_z,\
    int *flag_out,unsigned long long *N,unsigned long long *stepLineLen);

__global__ void TestMem(int *flag_out);