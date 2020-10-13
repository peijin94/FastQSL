import numpy as np
import scipy.io as sciIO

def trilerp(Bx,By,Bz,vec):
    idx0 = np.int32(np.floor(vec))
    w0 = np.stack(((1-(vec-idx0)),((vec-idx0))))
    weight = np.array([w0[i,0]*w0[j,1]*w0[k,2] for i in range(2) for j in range(2) for k in range(2)]).reshape(2,2,2) 
    bx = np.sum(Bx[idx0[0]:(idx0[0]+2),idx0[1]:(idx0[1]+2),idx0[2]:(idx0[2]+2)]*weight)
    by = np.sum(By[idx0[0]:(idx0[0]+2),idx0[1]:(idx0[1]+2),idx0[2]:(idx0[2]+2)]*weight)
    bz = np.sum(Bz[idx0[0]:(idx0[0]+2),idx0[1]:(idx0[1]+2),idx0[2]:(idx0[2]+2)]*weight)
    Bvec = np.array([bx,by,bz])
    return Bvec

