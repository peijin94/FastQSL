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


a = sciIO.readsav('../B0.sav')
print(a.keys())
Bx = a.twbox[0].bx.astype(np.float32)
By = a.twbox[0].by.astype(np.float32)
Bz = a.twbox[0].bz.astype(np.float32)
Bx = gpuarray.to_gpu(Bx)
By = gpuarray.to_gpu(By)
Bz = gpuarray.to_gpu(Bz)



# shape of B
BshapeN = np.zeros(3,dtype=np.int32)
BshapeN[:] = Bx.shape


interp_ratio=4
x_range = [200,400]
y_range = [200,400]
x_i = np.linspace(*x_range, interp_ratio*(x_range[1]-x_range[0]))
y_i = np.linspace(*y_range, interp_ratio*(x_range[1]-x_range[0]))

x_arr,y_arr = xx, yy = np.meshgrid(x_i, y_i)

xy_shape = x_arr.shape

x_inp = x_arr.flatten().astype(np.float32)
y_inp = y_arr.flatten().astype(np.float32)
z_inp = np.zeros_like(x_inp).astype(np.float32)

x_out = np.zeros_like(x_inp).astype(np.float32)
y_out = np.zeros_like(x_inp).astype(np.float32)
z_out = np.zeros_like(x_inp).astype(np.float32)


flag_out = np.zeros_like(x_inp).astype(np.int32)

s_len = np.float32([0.25])
N=np.ulonglong([x_inp.shape[0]])

blck = (16,16,1)
grid_a = 2**(np.int(np.log2(np.sqrt(N/256))))
grd = (grid_a,int(np.ceil(N/grid_a/256)))


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

pycuda.driver.Context.synchronize()
TraceAllBline(Bx,By,Bz,BshapeN,
             x_inp,y_inp,z_inp,
             x_out,y_out,z_out,
             s_len,flag_out,N,
             block=blck,grid=grd)
pycuda.driver.Context.synchronize()

print('haha')
print(blck)
print(grd)

print(N.get())
print(x_out.get())