import numpy as np
import cupy

def CookPseudoLine(x_end_arr,   y_end_arr,  z_end_arr,   flag_end_arr,
     x_start_arr, y_start_arr,z_start_arr, flag_start_arr,
     Bx_in_arr,   By_in_arr,  Bz_in_arr,
     Bx_out_arr,  By_out_arr, Bz_out_arr,
     Bz_0_arr,    B_flag_arr, stride_this):
    
    # collect index for pseudo lines
    (Bx_cut_start,By_cut_start,Bz_cut_start,flag_inp_cut,
    Bz_cut_end,Bz_cut_this,x_inp_cut0,y_inp_cut0,z_inp_cut0,
    x_out_cut0,y_out_cut0,z_out_cut0) = [arr[B_flag_arr==1] for arr in (
    Bx_in_arr,By_in_arr,Bz_in_arr,flag_start_arr,
    Bz_out_arr,Bz_0_arr,x_start_arr,y_start_arr,z_start_arr,
    x_end_arr,y_end_arr,z_end_arr)]
    
    Bz0_start=Bz_cut_start; Bz0_end=Bz_cut_end;
    
    # init gpu mem
    (cut_inp_x,cut_inp_y,cut_inp_z,cut_start_x,cut_start_y,cut_start_z,
     cut_end_x,cut_end_y,cut_end_z,Bx_start_cut,By_start_cut,Bz_start_cut,
     Bx_end_cut,By_end_cut,Bz_end_cut,Bx_inp_cut,By_inp_cut,Bz_inp_cut
    ) = [cupy.zeros(Bx_cut_start.shape[0]*4,dtype=cupy.float32) for _ in range(18)]
    line_len_cut = cupy.zeros(Bx_cut_start.shape[0]*4,cupy.float64)
    N_cut = cupy.array([cut_inp_z.shape[0]],cupy.ulonglong)
    (B_flag_cut,flag_cut_start,flag_cut_end)=[
        cupy.zeros(cut_inp_x.shape,dtype=cupy.int32) for _ in range(3)]
    
    # indexing
    (idx_X_cut,idx_Y_cut,idx_Z_cut) = [((flag_inp_cut-1)//2 == i) for i in range(3)]
    for arr_prep,arr_src in zip([cut_inp_x,cut_inp_y,cut_inp_z],
                    [x_inp_cut0,y_inp_cut0,z_inp_cut0]):
        for idx_start in range(4):
            arr_prep[idx_start::4]=arr_src
    
    # Z Y X direction
    cut_inp_x[cupy.where(idx_Z_cut)[0]*4+0]=x_inp_cut0[idx_Z_cut]+stride_this
    cut_inp_x[cupy.where(idx_Z_cut)[0]*4+2]=x_inp_cut0[idx_Z_cut]-stride_this
    cut_inp_y[cupy.where(idx_Z_cut)[0]*4+1]=y_inp_cut0[idx_Z_cut]+stride_this
    cut_inp_y[cupy.where(idx_Z_cut)[0]*4+3]=y_inp_cut0[idx_Z_cut]-stride_this
    cut_inp_x[cupy.where(idx_Y_cut)[0]*4+0]=x_inp_cut0[idx_Y_cut]+stride_this
    cut_inp_x[cupy.where(idx_Y_cut)[0]*4+2]=x_inp_cut0[idx_Y_cut]-stride_this
    cut_inp_z[cupy.where(idx_Y_cut)[0]*4+1]=z_inp_cut0[idx_Y_cut]+stride_this
    cut_inp_z[cupy.where(idx_Y_cut)[0]*4+3]=z_inp_cut0[idx_Y_cut]-stride_this
    cut_inp_y[cupy.where(idx_X_cut)[0]*4+0]=y_inp_cut0[idx_X_cut]+stride_this
    cut_inp_y[cupy.where(idx_X_cut)[0]*4+2]=y_inp_cut0[idx_X_cut]-stride_this
    cut_inp_z[cupy.where(idx_X_cut)[0]*4+1]=z_inp_cut0[idx_X_cut]+stride_this
    cut_inp_z[cupy.where(idx_X_cut)[0]*4+3]=z_inp_cut0[idx_X_cut]-stride_this

    return (cut_inp_x,   cut_inp_y,   cut_inp_z,
            cut_start_x, cut_start_y, cut_start_z,flag_cut_start,
            cut_end_x,   cut_end_y,   cut_end_z,  flag_cut_end,
            Bx_inp_cut,  By_inp_cut,  Bz_inp_cut,  B_flag_cut,
            Bx_start_cut,By_start_cut,Bz_start_cut,
            Bx_end_cut,  By_end_cut,  Bz_end_cut,  N_cut,line_len_cut,
            Bz0_start, Bz0_end)

def ResReshape(shape,*args):
    res_list=[]
    for res_var in args:
        res_list.append(res_var.reshape(shape))
    return res_list

def QCalcPlane(x_end_arr,   y_end_arr,  z_end_arr,   flag_end_arr,
     x_start_arr, y_start_arr,z_start_arr, flag_start_arr,
     Bx_in_arr,   By_in_arr,  Bz_in_arr,
     Bx_out_arr,  By_out_arr, Bz_out_arr,
     Bz_0_arr,    B_flag_arr, stride_step):

    (X1,Y1,X2,Y2)= [cupy.zeros(B_flag_arr.shape,dtype=cupy.float32) 
                   for _ in range(4)];    
    
    (idx_X1,idx_Y1,idx_Z1) = [((flag_start_arr-1)//2 == i) for i in range(3)]
    (idx_X2,idx_Y2,idx_Z2) = [( (flag_end_arr-1)//2  == i) for i in range(3)]

    X1[idx_Z1] = x_start_arr[idx_Z1]; Y1[idx_Z1] = y_start_arr[idx_Z1] 
    X2[idx_Z2] = x_end_arr[idx_Z2];   Y2[idx_Z2] = y_end_arr[idx_Z2] 
    X1[idx_Y1] = z_start_arr[idx_Y1]; Y1[idx_Y1] = x_start_arr[idx_Y1]
    X2[idx_Y2] = z_end_arr[idx_Y2];   Y2[idx_Y2] = x_end_arr[idx_Y2]
    X1[idx_X1] = y_start_arr[idx_X1]; Y1[idx_X1] = z_start_arr[idx_X1]
    X2[idx_X2] = y_end_arr[idx_X2];   Y2[idx_X2] = z_end_arr[idx_X2]

    dx2xc = X2[2:,1:-1]-X2[0:-2,1:-1]; dx2yc = X2[1:-1,2:]-X2[1:-1,0:-2];
    dy2xc = Y2[2:,1:-1]-Y2[0:-2,1:-1]; dy2yc = Y2[1:-1,2:]-Y2[1:-1,0:-2];
    dx1xc = X1[2:,1:-1]-X1[0:-2,1:-1]; dx1yc = X1[1:-1,2:]-X1[1:-1,0:-2];
    dy1xc = Y1[2:,1:-1]-Y1[0:-2,1:-1]; dy1yc = Y1[1:-1,2:]-Y1[1:-1,0:-2];

    a = (dx2xc*dy1yc-dx2yc*dy1xc);
    b = (dx2yc*dx1xc-dx2xc*dx1yc);
    c = (dy2xc*dy1yc-dy2yc*dy1xc);
    d = (dy2yc*dx1xc-dy2xc*dx1yc);

    bnr = cupy.abs(Bz_in_arr[1:-1,1:-1])/(cupy.abs(Bz_0_arr[1:-1,1:-1])**2
                    ) *cupy.abs(Bz_out_arr[1:-1,1:-1]) *((1/stride_step/2)**4)
    Q = (a**2+b**2+c**2+d**2)*bnr
    Q[cupy.where(Q<1.0)]=1.0
    return Q


def trilerp(Bx,By,Bz,vec):
    idx0 = np.int32(np.floor(vec))
    w0 = np.stack(((1-(vec-idx0)),((vec-idx0))))
    weight = np.array([w0[i,0]*w0[j,1]*w0[k,2] for i in range(2) 
                       for j in range(2) for k in range(2)]).reshape(2,2,2) 
    bx = np.sum(Bx[idx0[0]:(idx0[0]+2),idx0[1]:(idx0[1]+2),idx0[2]:(idx0[2]+2)]*weight)
    by = np.sum(By[idx0[0]:(idx0[0]+2),idx0[1]:(idx0[1]+2),idx0[2]:(idx0[2]+2)]*weight)
    bz = np.sum(Bz[idx0[0]:(idx0[0]+2),idx0[1]:(idx0[1]+2),idx0[2]:(idx0[2]+2)]*weight)
    Bvec = np.array([bx,by,bz])
    return Bvec

