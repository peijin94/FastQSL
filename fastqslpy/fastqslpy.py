import os

# get the path of the current file
path = os.path.dirname(os.path.realpath(__file__))
src_cu_path = os.path.join(path, 'src', '/cu/')


def kernel_TraceBlineAdaptive():
    print('compiling kernel')
    traceFunc_file = open(src_cu_path+"TraceBlineAdaptive.cu", "rt")
    TraceFunc =cupy.RawModule(code=traceFunc_file.read(),backend='nvcc',options=("-I "+src_cu_path,))#, include_dirs=[PWD], cache_dir='cache',no_extern_c=True)
    TraceAllBline = TraceFunc.get_function("TraceAllBline")
    return TraceAllBline

    