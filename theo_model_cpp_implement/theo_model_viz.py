
import numpy as np
import matplotlib.pyplot as plt
theomodel_path = "./build/rtm_theo_model.txt"
theoretical_model = np.loadtxt(theomodel_path)

cx = np.array((-1,1))
Lx,Ly = 2000,1200
plt.imshow(theoretical_model.T, aspect='auto',cmap='jet',extent=[0,Lx,Ly,0])
plt.clim(cx * 5e-16)
plt.colorbar()
plt.title("2-layer heterogeneous model")
plt.show()
