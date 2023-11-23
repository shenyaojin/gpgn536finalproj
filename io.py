import numpy as np
model= np.load('marmousi.npy')
# I/O to convert encrypted .npy file to .txt file
file_path = './marmousi_cppimplement/data.txt'
# Save the 2D array to a text file
model = model.T
np.savetxt(file_path, model, fmt='%.4f', delimiter='\t')