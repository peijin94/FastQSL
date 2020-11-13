import numpy as np
import pyvista as pv
npzfile = np.load('../QubePOT.npz')
Qube = npzfile['Qube']
#Twube = npzfile['arr_4']

axis_r = 0.364257
xi = npzfile['xi']*axis_r
yi = npzfile['yi']*axis_r
zi = npzfile['zi']*axis_r

Bz0 = npzfile['Bz0']


xx3,yy3,zz3 = np.meshgrid(xi, yi, zi,indexing='ij')

xx2,yy2 = np.meshgrid(xi, yi,indexing='ij')
zz2 = np.zeros_like(xx2)
B0mesh = pv.StructuredGrid(xx2, yy2, zz2)

grid = pv.StructuredGrid(xx3, yy3, zz3)
grid["vol"] = Qube.T.flatten()
contours = grid.contour([40000])

#grid["vol"] = Twube.T.flatten()
#contours = grid.contour([-10])

clim = [4000, 40000]
p = pv.Plotter()
largest = contours.connectivity(largest=True)
surf = largest.extract_geometry()
#surf = contours.extract_geometry()

#smooth = surf.smooth(n_iter=10)
#pv.set_plot_theme('document')
p = pv.Plotter()
#9ac4d9
p.add_mesh(surf, color='#9ac4d9', smooth_shading=True,opacity=0.4 )
p.add_mesh(B0mesh,scalars=Bz0.T, cmap='gray',clim=[-1500,1500],show_scalar_bar=True)
#p.enable_eye_dome_lighting()
#p.add_mesh(surf.outline(), color="k")
p.camera_position = [(309.21929818544663, -39.476961535608346, 62.5609895585235),
 (206.8292875784721, 95.43630460720914, 14.52431785900496),
 (-0.16273952375365774, 0.21903523353563947, 0.9620495901347398)]
p.show_bounds()
a=p.show()


print(a)

