import numpy as np
import pyvista as pv
npzfile = np.load('../Qube.npz')
Qube = npzfile['arr_0']
xx = npzfile['arr_1']
yy = npzfile['arr_2']
zz = npzfile['arr_3']

grid = pv.StructuredGrid(xx, yy, zz)

grid["vol"] = Qube.T.flatten()
contours = grid.contour([4000])
clim = [4000, 40000]
p = pv.Plotter()
largest = contours.connectivity(largest=True)
surf = largest.extract_geometry()
#smooth = surf.smooth(n_iter=10)
#pv.set_plot_theme('document')
p = pv.Plotter()
#9ac4d9
p.add_mesh(surf, color='#9ac4d9', smooth_shading=True,opacity=0.4 )
#p.enable_eye_dome_lighting()
p.add_mesh(surf.outline(), color="k")
p.show_bounds()
p.show()

