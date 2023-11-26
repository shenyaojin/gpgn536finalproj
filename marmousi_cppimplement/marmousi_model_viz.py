
import numpy as np
import matplotlib.pyplot as plt
marmousimodel_path = "./build/rtm_marmousi_model.txt"
marmousimodel = np.loadtxt(marmousimodel_path)

# plot parameters for marmousi model
nz,nx = 300,298
dz,dx = 5.0,5.0
oz,ox = 0.0,0.0
z = np.linspace(0,(nz-1)*dz,nz)
x = np.linspace(0,(nx-1)*dx,nx)
Lz, Lx = nz * dz, nx * dx
cx = np.array((-1,1))
plt.figure()
plt.imshow(marmousimodel.T, aspect='auto',cmap='jet',extent=[0,Lx,Lz,0])
plt.clim(cx * 1e-15)
plt.title("RTM result of Marmousi Model(nt = 2500)")
plt.xlabel("X/m")
plt.ylabel("Depth/m")
plt.colorbar()
plt.show()
