import numpy as np
import pyvista as pv
npzfile = np.load('../QubePOT.npz')
Qube = npzfile['Qube']
#Twube = npzfile['arr_4']

xi = npzfile['xi']
yi = npzfile['yi']
zi = npzfile['zi']

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
p.camera_position = [(918.4679349260575, -222.2585068804345, 282.83224256932135),
 (510.0, 320.0, 60.0),
 (-0.13889900119010723, 0.28509344203728243, 0.9483821997358055)]
p.show_bounds()
a=p.show()


print(a)

