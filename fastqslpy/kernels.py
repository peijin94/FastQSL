import os
import cupy

# get the path of the current file
dirname = os.path.dirname(os.path.abspath(__file__))
src_cu_path = os.path.join(dirname, 'cu')

modules = ['TraceBline','TraceBlineAdaptive','TraceBlineAdaptiveStretch','TraceBlineScott']

def compileTraceBlineAdaptive():
    print('compiling kernel')
    print(dirname)
    print(__file__)
    print(src_cu_path)
    traceFunc_file = open(os.path.join(src_cu_path,"TraceBlineAdaptive.cu"), "rt")
    TraceFunc =cupy.RawModule(code=traceFunc_file.read(),backend='nvcc',options=("-I "+src_cu_path,))#, include_dirs=[PWD], cache_dir='cache',no_extern_c=True)
    TraceAllBline = TraceFunc.get_function("TraceAllBline")
    return TraceAllBline
