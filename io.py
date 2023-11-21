import numpy as np
model= np.load('marmousi.npy')
# I/O to convert encrypted .npy file to .txt file
file_path = './cppimplement/data.txt'
# Save the 2D array to a text file
np.savetxt(file_path, model, fmt='%.4f', delimiter='\t')