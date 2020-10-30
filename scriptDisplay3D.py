import numpy as np
import pyvista as pv
npzfile = np.load('../Qube.npz')
Qube = npzfile['arr_0']
xx = npzfile['arr_1']
yy = npzfile['arr_2']
zz = npzfile['arr_3']

grid = pv.StructuredGrid(xx.flatten(), yy.flatten(), zz.flatten())

grid["vol"] = Qube.flatten()
contours = grid.contour([100])
clim = [4000, 40000]
p = pv.Plotter()
largest = contours.connectivity(largest=True)
pv.set_plot_theme('document')
p = pv.Plotter()

p.add_mesh(contours, color='#965434')
p.enable_eye_dome_lighting()
#p.add_mesh(mesh.outline(), color="k")
#p.add_mesh(largest, scalars=contours.points[:, 2], show_scalar_bar=False)
p.show()


#p.add_volume(Qube, cmap="inferno", clim=clim,
#             opacity=opacity, opacity_unit_distance=4,mapper='gpu')
#p.camera_position =[(-800, 1200, 666),
# (179.5, 299.5, 99.5),
# (0.4, -0.1, 0.9)]
#p.add_bounds_axes()
#p.show_axes()
#a=p.show()
#print(p)
#print(a)