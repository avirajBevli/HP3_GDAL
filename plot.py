import os
import subprocess
from osgeo import gdal
import matplotlib.pyplot as plt

slp1 = gdal.Open("cdnf43w.tif")
slp1Array = slp1.GetRasterBand(1).ReadAsArray()

plt.figure()
plt.imshow(slp1Array)
plt.colorbar()
plt.show()