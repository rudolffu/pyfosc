#!/usr/bin/env python
# read in the standard star spectrum with the following format from dat file:
# wavelength AB_mag binwidth(optional)
# 1147.9000 12.882  1.2
# 1149.1000 12.988  1.2
# 1150.3000 12.915  1.2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
import os
import sys

filename = sys.argv[1]
abspath = os.path.abspath(filename)
output_path = os.path.dirname(os.path.dirname(abspath))
objname = sys.argv[2]
# check if the file has 3 columns
with open(filename, 'r') as f:
    first_line = f.readline()
    
if len(first_line.split()) == 2:
    df = pd.read_csv(filename, sep='\s+',
                        header=None,
                        names=['wave', 'AB_mag'])
elif len(first_line.split()) == 3:
    df = pd.read_csv(filename, sep='\s+',
                        header=None,
                        names=['wave', 'AB_mag', 'binwidth'])

# resample the spectrum to 20-Angstrom bins
from scipy.interpolate import interp1d
f = interp1d(df.wave, df.AB_mag)
xnew = np.arange(df.wave.min(), df.wave.max(), 20)
ynew = f(xnew)
# round all values to 3 decimal places
xnew = np.round(xnew, 3)
ynew = np.round(ynew, 3)
xnew = np.concatenate((xnew, [df.wave.values[-1]]))
ynew = np.concatenate((ynew, [df.AB_mag.values[-1]]))
wnew = 20 * np.ones_like(xnew)

# plot the resampled spectrum
plt.figure()
plt.plot(xnew, ynew)
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('AB mag')
plt.title('Resampled Spectrum')
plt.show()

# save the resampled spectrum to a file
df_new = pd.DataFrame({'wave': xnew, 'AB_mag': ynew, 'binwidth': wnew})
output_filename = os.path.join(output_path, objname + '.dat')
df_new.to_csv(output_filename, sep='\t', index=False, header=False)