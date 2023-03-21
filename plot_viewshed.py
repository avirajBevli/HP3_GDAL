import matplotlib.pyplot as plt
import numpy as np
  
# Text file data converted to integer data type
File_data = np.loadtxt("viewshed_r3.txt", dtype=int)
# print(File_data)
plt.figure()
plt.imshow(File_data)
plt.colorbar()
plt.show()