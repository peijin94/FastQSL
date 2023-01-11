import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import pycuda.gpuarray as gpuarray
import numpy as np
import pycuda.driver as cuda
import scipy.io as sciIO


# compile from CUDA/C++ code
traceFunc_file = open("TraceBline.cu", "rt")
traceFunc =SourceModule(traceFunc_file.read()) 

TraceAllBline = traceFunc.get_function("TraceAllBline")
test_Idx3d = traceFunc.get_function("test_Idx3d")
test_Interp3d = traceFunc.get_function("test_Interp3d")


a = sciIO.readsav('../B0.sav')
print(a.keys())
Bx = np.transpose(a.twbox[0].bx.astype(np.float32),(2,1,0))
By = np.transpose(a.twbox[0].by.astype(np.float32),(2,1,0))
Bz = np.transpose(a.twbox[0].bz.astype(np.float32),(2,1,0))
Bx_gpu = gpuarray.to_gpu(Bx)
By_gpu = gpuarray.to_gpu(By)
Bz_gpu = gpuarray.to_gpu(Bz)



realpos=[586,300,26]
print("Host Value : " +str(Bz[realpos[0],realpos[1],realpos[2]]))

# shape of B
BshapeN = np.zeros(3,dtype=np.int_)
BshapeN[:] = Bz.shape

#getIdx = (1.1+np.array(realpos)).astype(np.float32)
getIdx = np.array(realpos).astype(np.int_)
res = np.array([0.0]).astype(np.float32)

print('Index : '+str(getIdx))
print('Shape : '+str(BshapeN))

Bz_gpu = gpuarray.to_gpu(Bz)
getIdx_gpu = gpuarray.to_gpu(getIdx)
res_gpu = gpuarray.to_gpu(res)
BshapeN = gpuarray.to_gpu(BshapeN)

# test indexing 
test_Idx3d(Bz_gpu,BshapeN,getIdx_gpu,res_gpu,block=(1,1,1))
print("GPU  Value : " +str(res_gpu.get()))
cuda.Context.synchronize()

# test interpolation
pinPoint = (np.array([586.97,300.61,25.58])).astype(np.float32)
pinPoint_gpu = gpuarray.to_gpu(pinPoint)
print(pinPoint)
test_Interp3d(Bz_gpu,BshapeN,pinPoint_gpu,res_gpu,block=(1,1,1))
print(res_gpu.get())
cuda.Context.synchronize()
