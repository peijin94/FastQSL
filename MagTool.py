import numpy as np
import cupy

# use cupy as numker if Nvidia GPU exsit
from numba import jit

#@jit(fastmath = True)
def trilerp(Bx,By,Bz,vec,numker=np):
    vec[vec<0]=0
    idx=0
    for val in vec:
        if val>Bx.shape[idx]:
            vec[idx]=Bx.shape[idx]-1e-8
        idx=idx+1
    idx0 = numker.floor(vec).astype(numker.int)
    Bxblock=Bx[idx0[0]:(idx0[0]+2),idx0[1]:(idx0[1]+2),idx0[2]:(idx0[2]+2)]
    Byblock=By[idx0[0]:(idx0[0]+2),idx0[1]:(idx0[1]+2),idx0[2]:(idx0[2]+2)]
    Bzblock=Bz[idx0[0]:(idx0[0]+2),idx0[1]:(idx0[1]+2),idx0[2]:(idx0[2]+2)]
    w0 = numker.stack(((1-(vec-idx0)),((vec-idx0))))
    weight= numker.array([w0[i,0]*w0[j,1]*w0[k,2] for i in range(2) 
                       for j in range(2) for k in range(2)]).reshape(2,2,2)
    bx = numker.sum(Bxblock*weight)
    by = numker.sum(Byblock*weight)
    bz = numker.sum(Bzblock*weight)
    Bvec = numker.array([bx,by,bz])
    return Bvec


def TraceBline(Bx,By,Bz,Pxyz0,dr=0.25,maxstep=1e4,numker=np):
    maxstep=numker.int(maxstep)
    B0 = trilerp(Bx,By,Bz,Pxyz0,numker=numker)
    Direction  = numker.sign(B0[2])
    def Bnorm(Bx,By,Bz,vec,Direction=Direction,numker=numker):
        B0=trilerp(Bx,By,Bz,vec,numker=numker)
        return B0/numker.sqrt(numker.sum(B0**2))*Direction
    PP0 = Pxyz0
    (xline,yline,zline)=[numker.zeros(maxstep,dtype=numker.float32)+numker.nan for _ in range(3)]
    xline[0] = Pxyz0[0]
    yline[0] = Pxyz0[1]
    zline[0] = Pxyz0[2]
    # use RK4 to solve r
    for idx in range(maxstep):
        f1 = Bnorm(Bx,By,Bz,PP0)
        f2 = Bnorm(Bx,By,Bz,PP0+dr*f1/2.)
        f3 = Bnorm(Bx,By,Bz,PP0+dr*f2/2.)
        f4 = Bnorm(Bx,By,Bz,PP0+dr*f3)
        xline[idx] = PP0[0]
        yline[idx] = PP0[1]
        zline[idx] = PP0[2]
        PP0 = PP0+dr*(f1+2.*f2+2.*f3+f4)/6.
        if ( PP0[0]<0 or PP0[1]<0 or PP0[2]<0 or 
            PP0[0]>Bx.shape[0] or PP0[1]>Bx.shape[0] or PP0[2]>Bx.shape[0]):
            break
    return (xline,yline,zline)
    
    