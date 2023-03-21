import os
import subprocess
from osgeo import gdal
import matplotlib.pyplot as plt

slp1 = gdal.Open("cdnf43w.tif")
slp1Array = slp1.GetRasterBand(1).ReadAsArray()

print(type(slp1Array))
print(slp1Array.shape) # (3600,3600)


# open a file for writing and write text in it
with open('heights.txt', 'w') as f:
    for i in range(slp1Array.shape[0]):
    	for j in range(slp1Array.shape[1]):
    		f.write(str(slp1Array[i][j])) # need to convert into string before writing to file
    		f.write(', ')
    	f.write('\n')

plt.figure()
plt.imshow(slp1Array)
plt.colorbar()
plt.show()