import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import pycuda.gpuarray as gpuarray
import numpy as np
import pycuda.driver as cuda

import scipy.io as sciIO

traceFunc_file = open("TraceBline.cu", "rt")
traceFunc =SourceModule(traceFunc_file.read())

Interp3d = traceFunc.get_function("Interp3d")

a = sciIO.readsav('../B0.sav')
print(a.keys())
Bx = a.twbox[0].bx.astype(np.float32)
By = a.twbox[0].by.astype(np.float32)
Bz = a.twbox[0].bz.astype(np.float32)
#Bx = gpuarray.to_gpu(Bx)
#By = gpuarray.to_gpu(By)
Bz = gpuarray.to_gpu(Bz)

BshapeN = np.zeros(3,dtype=np.int_)
BshapeN[:] = Bx.shape
getIdx = np.zeros(3,dtype=np.int_)
getIdx[:] = np.array([9,2,1],dtype=np.int_)
res=np.array([0.0],dtype=np.float32)
getPoint = np.array([1.1,1.2,1.1],dtype=np.float32)
#garray.empty((1), dtype=np.float)


BshapeN = gpuarray.to_gpu(BshapeN)
getIdx = gpuarray.to_gpu(getIdx)
res_gpu = gpuarray.to_gpu(res)
getPoint_gpu = gpuarray.to_gpu(getPoint)

Interp3d(Bz,BshapeN,getPoint_gpu,res_gpu,block=(1,1,1))

print(res_gpu.get())