import numpy as np

def trilerp(Bx,By,Bz,vec,nk=np): #interpolate B_vec from Bx,By,Bz
    idx=0
    for val in vec:
        if val<0:
            vec[idx]=1e-8
        if val>Bx.shape[idx]:
            vec[idx]=Bx.shape[idx]-1e-8
        idx=idx+1
    idx0 = (vec).astype(nk.int)
    Bxblock=Bx[idx0[0]:(idx0[0]+2),idx0[1]:(idx0[1]+2),idx0[2]:(idx0[2]+2)]
    Byblock=By[idx0[0]:(idx0[0]+2),idx0[1]:(idx0[1]+2),idx0[2]:(idx0[2]+2)]
    Bzblock=Bz[idx0[0]:(idx0[0]+2),idx0[1]:(idx0[1]+2),idx0[2]:(idx0[2]+2)]
    w0 = nk.stack(((1-(vec-idx0)),((vec-idx0))))
    weight= nk.array([w0[i,0]*w0[j,1]*w0[k,2] for i in range(2) 
                       for j in range(2) for k in range(2)]).reshape(2,2,2)
    return nk.array([nk.sum(B*weight) for B in (Bxblock,Byblock,Bzblock)])
    
def Bnorm(Bx,By,Bz,vec,nk=np):
    B0=trilerp(Bx,By,Bz,vec,nk=nk)
    return B0/nk.sqrt(nk.sum(B0**2))

def RK4_bline(Bx,By,Bz,PP0,xl,yl,zl,steprange,dr=0.25,nk=np):
    for idx in steprange:
        if ( PP0[0]<0          or PP0[1]<0           or PP0[2]<0 or 
            PP0[0]>Bx.shape[0] or PP0[1]>Bx.shape[0] or PP0[2]>Bx.shape[0]):
            break
        f1 = Bnorm(Bx,By,Bz,PP0)
        f2 = Bnorm(Bx,By,Bz,PP0+dr*f1/2.)
        f3 = Bnorm(Bx,By,Bz,PP0+dr*f2/2.)
        f4 = Bnorm(Bx,By,Bz,PP0+dr*f3)
        (xl[idx],yl[idx],zl[idx])=PP0
        PP0 = PP0+dr*(f1+2.*f2+2.*f3+f4)/6.

def TraceBline(Bx,By,Bz,Pxyz0,dr=0.25,maxstep=2e4,nk=np):
    maxstep=nk.int(maxstep)
    midpoint=nk.int(maxstep/2) 
    (xline,yline,zline)=[nk.zeros(maxstep,dtype=nk.float32)+nk.nan for _ in range(3)]
    # use RK4 to solve r
    PP0 = Pxyz0
    RK4_bline(Bx,By,Bz,PP0,xline,yline,zline,range(midpoint,maxstep),dr=0.25,nk=nk) # forward
    RK4_bline(Bx,By,Bz,PP0,xline,yline,zline,range(midpoint,0,-1),  dr=-0.25,nk=nk) # backward
    return (xline,yline,zline)