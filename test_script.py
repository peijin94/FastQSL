import numpy as np
import scipy.io as sciIO
import matplotlib.pyplot as plt
import pycuda.gpuarray as gpuarray
import pycuda.driver as cuda

import time


a = sciIO.readsav('../B0.sav')
print(a.keys())

# need to be transposed because of the storage toplogy
Bx = np.transpose(a.twbox[0].bx,(2,1,0))
By = np.transpose(a.twbox[0].by,(2,1,0))
Bz = np.transpose(a.twbox[0].bz,(2,1,0))

Bx_gpu = np.zeros(0)
By_gpu = np.zeros(0)
Bz_gpu = np.zeros(0)

import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import numpy

traceFunc_file = open("TraceBline.cu", "rt")
traceFunc =SourceModule(traceFunc_file.read())
print('compiling kernel')
TraceAllBline = traceFunc.get_function("TraceAllBline")

if Bx_gpu.shape[0]<=1:
    print('transfering B-field to GPU')
    Bx_gpu = gpuarray.to_gpu(Bx.astype(np.float32))
    By_gpu = gpuarray.to_gpu(By.astype(np.float32))
    Bz_gpu = gpuarray.to_gpu(Bz.astype(np.float32))

# shape of B
BshapeN = np.zeros(3,dtype=np.int32)
BshapeN[:] = Bx.shape

interp_ratio=4
x_range = [0,1163]
y_range = [0,487]
#x_range = [350,670]
#y_range = [200,400]
#x_range = [450,550]
#y_range = [280,350]
#x_range = [587,588]
#y_range = [301,302]
x_i = np.linspace(*x_range, np.uint32(interp_ratio*(x_range[1]-x_range[0])))
y_i = np.linspace(*y_range, np.uint32(interp_ratio*(y_range[1]-y_range[0])))
x_arr,y_arr = np.meshgrid(x_i, y_i)

xy_shape = x_arr.shape

x_inp = x_arr.flatten().astype(np.float32)
y_inp = y_arr.flatten().astype(np.float32)
z_inp = np.zeros_like(x_inp).astype(np.float32)

#z_inp[:] = 23.0

x_out = np.zeros_like(x_inp).astype(np.float32)
y_out = np.zeros_like(x_inp).astype(np.float32)
z_out = np.zeros_like(x_inp).astype(np.float32)

Bx_out = np.zeros_like(x_inp).astype(np.float32)
By_out = np.zeros_like(x_inp).astype(np.float32)
Bz_out = np.zeros_like(x_inp).astype(np.float32)

Bx_inp = np.zeros_like(x_inp).astype(np.float32)
By_inp = np.zeros_like(x_inp).astype(np.float32)
Bz_inp = np.zeros_like(x_inp).astype(np.float32)

flag_out = np.zeros_like(x_inp).astype(np.int32)
step_line_len = np.zeros_like(x_inp).astype(np.ulonglong)

s_len = np.float32([0.25])
N=np.ulonglong([x_inp.shape[0]])

# for GTX 1060
#blck=(64,1,1)
#grd = (20,1)

# for GTX1080ti
#blck=(128,1,1)
#grd = (28,1)


# for GTX1080ti
blck=(128,1,1)
grd = (60,1)

# chuck everything

BshapeN = gpuarray.to_gpu(BshapeN)
x_inp = gpuarray.to_gpu(x_inp)
y_inp = gpuarray.to_gpu(y_inp)
z_inp = gpuarray.to_gpu(z_inp)

x_out = gpuarray.to_gpu(x_out)
y_out = gpuarray.to_gpu(y_out)
z_out = gpuarray.to_gpu(z_out)

s_len    = gpuarray.to_gpu(s_len)
flag_out = gpuarray.to_gpu(flag_out)
N        = gpuarray.to_gpu(N)


Bx_out = gpuarray.to_gpu(Bx_out)
By_out = gpuarray.to_gpu(By_out)
Bz_out = gpuarray.to_gpu(Bz_out)

Bx_inp = gpuarray.to_gpu(Bx_inp)
By_inp = gpuarray.to_gpu(By_inp)
Bz_inp = gpuarray.to_gpu(Bz_inp)

step_line_len = gpuarray.to_gpu(step_line_len)


start = time.time()

TraceAllBline(Bx_gpu,By_gpu,Bz_gpu,BshapeN,
             x_inp,y_inp,z_inp,
             x_out,y_out,z_out,
             s_len,
             Bx_inp,By_inp,Bz_inp,
             Bx_out,By_out,Bz_out,
             flag_out,N,step_line_len,
             block=blck,grid=grd)



pycuda.driver.Context.synchronize()


end = time.time()
print(end - start)
