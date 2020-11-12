import numpy as np
import pyvista as pv
npzfile = np.load('../Qube.npz')
Qube = np.log10(npzfile['arr_0'])
Twube = npzfile['arr_4']
xx = npzfile['arr_1']
yy = npzfile['arr_2']
zz = npzfile['arr_3']

clim = np.log10(np.array([30000, 10000000]))
p = pv.Plotter()

opacity = [ 0,0.1,0.3, 0.45, 0.8]
p.add_volume(Qube, cmap="magma", clim=clim,
             opacity="sigmoid_5", mapper='gpu')#opacity_unit_distance=2,
#p.camera_position =[(-3074.26058295243, 1945.199902031422, 1494.7777365540755),
# (644.2846023046784, 730.3817379032646, 405.007411077699),
# (0.23837588336151042, -0.1332940950461866, 0.9619821320884487)]
p.show_grid()
p.show_axes()
a=p.show()
print(p)
print(a)