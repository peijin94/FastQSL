import numpy as np
import pyvista as pv
import time
npzfile = np.load('../QubeFinal.npz',allow_pickle=True)
Qube = (npzfile['Qube'])
Twube = npzfile['Twube']
xi = npzfile['xi']
yi = npzfile['yi']
zi = npzfile['zi']
blines=npzfile['blines']
B0z = npzfile['B0z']

xx3,yy3,zz3 = np.meshgrid(xi, yi, zi,indexing='ij')
xx2,yy2 = np.meshgrid(xi, yi,indexing='ij')

colors=["#ff5e41",*(["#f4ff41"]*8),"#58FF41",*(["#5189ff"]*3)]
radius=[0.6,*([0.3]*8),0.5 ,*([0.7]*3)]

p = pv.Plotter(off_screen=True)
p.open_movie("movie/test.mp4")
pv.set_plot_theme("ParaView")
idx=0
for bline in blines:
    line = pv.lines_from_points(bline)
    tube = line.tube(radius=radius[idx])
    p.add_mesh(tube,smooth_shading=True,color=colors[idx]) 
    idx=idx+1   

grid = pv.StructuredGrid(xx3, yy3, zz3,indexing='ij')
grid["vol"] = Qube.T.flatten()
contours = grid.contour([5000])
largest = contours.connectivity(largest=True)
surf = largest.extract_geometry()
p.add_mesh(surf,color='#9a04d9', smooth_shading=True, show_scalar_bar=False,opacity=0.07 )


xxi=np.arange(190,310)
yyi=np.arange(190,310)
xx2,yy2 = np.meshgrid(xxi,yyi,indexing='ij')
zz2 = np.zeros_like(xx2)
B0mesh = pv.StructuredGrid(xx2, yy2, zz2)
p.add_mesh(B0mesh,scalars=B0z[xx2,yy2].T, cmap='gray',clim=[-2,2])

theta=0.3*np.pi#t_id
p.camera_position =[(250+200*np.sin(theta), 200+80*np.cos(theta), 160),
     (250, 250, 30),
     (0, 0, 1)]

p.window_size = 1920, 1080
p.show_grid()
p.show_axes()

p.show(auto_close=False)

thetas  = (0.55+np.linspace(0,1,240))*np.pi
idx_this=0
for t_cur in thetas:
    theta=t_cur#0.3*np.pi#t_id
    p.camera_position =[(250+180*np.sin(theta), 250+180*np.cos(theta), 160),
         (250, 250, 30),
         (0, 0, 1)]
    p.write_frame() 
    #p.render()
    #p.screenshot('movie/img'+str(idx_this).rjust(3,'0')+'.png')
    #idx_this=idx_this+1

zs = 160-np.linspace(0,150,50)
for z_cur in zs:
    p.camera_position =[(250+180*np.sin(theta), 250+180*np.cos(theta), z_cur),
         (250, 250, 30),
         (0, 0, 1)]
    p.write_frame() 
    #p.render()
    #p.screenshot('movie/img'+str(idx_this).rjust(3,'0')+'.png')
    #idx_this=idx_this+1

    
ys = np.linspace(250+180*np.sin(theta),250,200)
for y_cur in ys:
    p.camera_position =[(y_cur,  250+180*np.cos(theta), 10),
         (250, 250, 30),
         (0, 0, 1)]
    p.write_frame() 
    #p.render()
    #p.screenshot('movie/img'+str(idx_this).rjust(3,'0')+'.png')
    #idx_this=idx_this+1
    
        
fs = np.linspace(30,9,100)
for f_cur in fs:
    p.camera_position =[(250,  250+180*np.cos(theta), 10),
         (250, 250, f_cur),
         (0, 0, 1)]
    p.write_frame() 
    #p.render()
    #p.screenshot('movie/img'+str(idx_this).rjust(3,'0')+'.png')
    #idx_this=idx_this+1


thetas  = (np.linspace(0,1,500))*np.pi
r=180*np.cos(theta)
for t_cur in thetas:
    theta2=t_cur
    p.camera_position =[(250+r*np.sin(theta2), 250+r*np.cos(theta2), 10),
         (250, 250, 9),
         (0, 0, 1)]
    p.write_frame() 
    
    
    
zs = np.linspace(10,250,400)
for z_cur in zs:
    p.camera_position =[(250+r*np.sin(theta2), 250+r*np.cos(theta2), z_cur),
         (250, 250, 9),
         (0, 0, 1)]
    p.write_frame() 
    
p.close()

print(p)
print(a)